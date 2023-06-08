#ifndef CIVET_SAMPLING_H
#define CIVET_SAMPLING_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <utils/rng.h>
#include <algorithm>

namespace civet {

struct Distribution1D {
	Distribution1D(const float* f, int n) :
			func(f, f + n), cdf(n + 1) {
		cdf[0] = 0;
		for (int i = 1; i < n + 1; i++) {
			cdf[i] = cdf[i - 1] + func[i - 1] / n;
		}

		func_int = cdf[n];
		if (func_int == 0) {
			for (int i = 1; i < n + 1; i++) {
				cdf[i] = float(i) / float(n);
			}
		} else {
			for (int i = 1; i < n + 1; i++) {
				cdf[i] / func_int;
			}
		}
	}

	int count() const { return func.size(); }

	float sampleContinuous(float u, float* pdf, int* off = nullptr) const {
		int offset = findInterval(cdf.size(), [&](int index) { return cdf[index] <= u; });
		if (off) {
			*off = offset;
		}
		float du = u - cdf[offset];
		if ((cdf[offset + 1] - cdf[offset]) > 0) {
			du /= (cdf[offset + 1] - cdf[offset]);
		}

		if (*pdf) {
			*pdf = func[offset] / func_int;
		}

		return (offset + du) / count();
	}

	float sampleDiscrete(float u, float* pdf, float* u_remapped = nullptr) {
		int offset = findInterval(cdf.size(), [&](int index) { return cdf[index] <= u; });
		if (pdf) {
			*pdf = func[offset] / (func_int * count());
		}
		if (u_remapped) {
			*u_remapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		}
		return offset;
	}

	float discretePDF(int index) const {
		return func[index] / (func_int * count());
	}

	std::vector<float> func, cdf;
	float func_int;
};

inline Point2f rejectionSampleDisk(RNG& rng) {
	Point2f p;
	do {
		p.x = 1 - 2 * rng.uniformFloat();
		p.y = 1 - 2 * rng.uniformFloat();
	} while (p.x * p.x + p.y * p.y > 1);
	return p;
}

inline Vector3f uniformSampleHemisphere(const Point2f& u) {
	float z = u[0];
	float r = std::sqrt(std::max(0.f, 1.f - z * z));
	float phi = 2 * Pi * u[1];
	return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline float uniformHemispherePdf() {
	return Inv2Pi;
}

inline Vector3f uniformSampleSphere(const Point2f& u) {
	float z = 2.f * u[0] - 1.f;
	float r = std::sqrt(std::max(0.f, 1.f - z * z));
	float phi = 2.f * Pi * u[1];
	return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

inline float uniformSpherePdf() {
	return Inv4Pi;
}

inline Point2f uniformSampleDisk(const Point2f& u) {
	float r = std::sqrt(u[0]);
	float theta = 2 * Pi * u[1];
	return Point2f(r * std::cos(theta), r * std::sin(theta));
}

inline Point2f concentricSampleDisk(const Point2f& u) {
	Point2f u_offset = 2.f * u - Vector2f(1, 1);
	if (u_offset.x == 0 && u_offset.y == 0) {
		return Point2f(0, 0);
	}
	float theta, r;
	if (std::abs(u_offset.x) > std::abs(u_offset.y)) {
		r = u_offset.x;
		theta = PiOver4 * (u_offset.y / u_offset.x);
	} else {
		r = u_offset.y;
		theta = PiOver2 - PiOver4 * (u_offset.x / u_offset.y);
	}
	return r * Point2f(std::cos(theta), std::sin(theta));
}

inline Vector3f uniformSampleCone(const Point2f& u, float cos_theta_max) {
	float cos_theta = (1.f - u[0]) + u[0] * cos_theta_max;
	float sin_theta = std::sqrt(1.f - cos_theta * cos_theta);
	float phi = u[1] * 2 * Pi;
	return Vector3f(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
}

inline Vector3f UniformSampleCone(const Point2f& u, float cos_theta_max, const Vector3f& x, const Vector3f& y, const Vector3f& z) {
	float cos_theta = lerp(u[0], cos_theta_max, 1.f);
	float sin_theta = std::sqrt(1.f - cos_theta * cos_theta);
	float phi = u[1] * 2 * Pi;
	return std::cos(phi) * sin_theta * x + std::sin(phi) * sin_theta * y + cos_theta * z;
}

inline float uniformConePdf(float cos_theta_max) {
	return 1 / (2 * Pi * (1 - cos_theta_max));
}

inline Point2f uniformSampleTriangle(const Point2f& u) {
	float su0 = std::sqrt(u[0]);
	return Point2f(1 - su0, u[1] * su0);
}

inline Vector3f cosineSampleHemisphere(const Point2f& u) {
	Point2f d = concentricSampleDisk(u);
	float z = std::sqrt(std::max(0.0f, 1.f - d.x * d.x - d.y * d.y));
	return Vector3f(d.x, d.y, z);
}

inline float cosineHemispherePdf(float cos_theta) {
	return cos_theta * InvPi;
}

inline float balanceHeuristic(int nf, float fpdf, int ng, float gpdf) {
	return (nf * fpdf) / (nf * fpdf + ng * gpdf);
}

inline float powerHeuristic(int nf, float fpdf, int ng, float gpdf) {
	float f = nf * fpdf, g = ng * gpdf;
	return (f * f) / (f * f + g * g);
}

inline float radicalInverseBase2(uint32_t bits) {
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	return float(bits) * 2.3283064365386963e-10f; // / 0x100000000
}

inline Vector2f hammersley2D(uint64_t sample_idx, uint64_t num_samples) {
	return Vector2f(float(sample_idx) / float(num_samples), radicalInverseBase2(uint32_t(sample_idx)));
}

// Generates hammersley using base 1 and 2
inline void generateHammersleySamples2D(Vector2f* samples, uint64_t numSamples) {
	for (uint64_t i = 0; i < numSamples; ++i) {
		samples[i] = hammersley2D(i, numSamples);
	}
}

} // namespace civet

#endif // CIVET_SAMPLING_H
