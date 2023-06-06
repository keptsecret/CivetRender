/*
 * Spherical gaussians and their usage adapted from: https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/SG.h
 * (MIT)
 */

#ifndef CIVET_SG_H
#define CIVET_SG_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>

namespace civet {

// SphericalGaussian(dir) := Amplitude * exp(Sharpness * (dot(Axis, Direction) - 1.0f))
struct SG {
	Vector3f amplitude;
	float sharpness = 1.0f;
	Vector3f axis;

	// exp(2 * Sharpness * (dot(Axis, Direction) - 1.0f)) integrated over the sampling domain.
	float basis_sq_integral_over_domain;
};

// Evaluates an SG given a direction on a unit sphere
inline Vector3f evaluteSG(const SG& sg, Vector3f dir) {
	return sg.amplitude * std::exp(sg.sharpness * (dot(sg.axis, dir) - 1.f));
}

// Computes the inner product of two SG's, which is equal to Integrate(SGx(v) * SGy(v) * dv).
inline Vector3f innerProduceSG(const SG& x, const SG& y) {
	float um_length = (x.sharpness * x.axis + y.sharpness * y.axis).length();
	Vector3f xaya(x.amplitude.x * y.amplitude.x, x.amplitude.y * y.amplitude.y, x.amplitude.z * y.amplitude.z);
	Vector3f expo = std::exp(um_length - x.sharpness - y.sharpness) * xaya;
	float other = 1.0f - std::exp(-2.0f * um_length);
	return (2.0f * Pi * expo * other) / um_length;
}

// Returns an approximation of the clamped cosine lobe represented as an SG
inline SG cosineLobeSG(Vector3f direction) {
	SG cosine_lobe;
	cosine_lobe.axis = direction;
	cosine_lobe.sharpness = 2.133f;
	cosine_lobe.amplitude = Vector3f(1.17f, 1.17f, 1.17f);
	return cosine_lobe;
}

// Computes the approximate integral of an SG over the entire sphere. The error vs. the
// non-approximate version decreases as sharpeness increases.
inline Vector3f approximateSGIntegral(const SG& sg) {
	return 2 * Pi * (sg.amplitude / sg.sharpness);
}

// Computes the approximate incident irradiance from a single SG lobe containing incoming radiance.
// The irradiance is computed using a fitted approximation polynomial. This approximation
// and its implementation were provided by Stephen Hill.
inline Vector3f SGIrradianceFitted(const SG& lighting_lobe, const Vector3f& normal) {
	const float muDotN = dot(lighting_lobe.axis, normal);
	const float lambda = lighting_lobe.sharpness;

	const float c0 = 0.36f;
	const float c1 = 1.0f / (4.0f * c0);

	float eml = std::exp(-lambda);
	float em2l = eml * eml;
	float rl = 1.0f / lambda;

	float scale = 1.0f + 2.0f * em2l - rl;
	float bias = (eml - em2l) * rl - em2l;

	float x = std::sqrt(1.0f - scale);
	float x0 = c0 * muDotN;
	float x1 = c1 * x;

	float n = x0 + x1;

	float y = (std::abs(x0) <= x1) ? n * n / x : clamp(muDotN, 0, 1);

	float normalized_irradiance = scale * y + bias;

	return normalized_irradiance * approximateSGIntegral(lighting_lobe);
}

// Input parameters for the solve
struct SGSolveParam {
	// StrikePlate plate;                           // radiance over the sphere
	Vector3f* x_samples = nullptr;
	Vector3f* y_samples = nullptr;
	uint64_t num_samples = 0;

	uint64_t num_SGs = 0; // number of SG's we want to solve for

	SG* out_SGs; // output of final SG's we solve for
};

enum class SGDistribution : uint32_t {
	Spherical,
	Hemispherical,
};

void initializeSGSolver(uint64_t numSGs, SGDistribution distribution);
const SG* initialGuess();

// Solve for k-number of SG's based on a hemisphere of radiance
void solveSGs(SGSolveParam& params);

void projectOntoSGs(const Vector3f& dir, const Vector3f& color, SG* out_SGs, uint64_t num_SGs);

void SGRunningAverage(const Vector3f& dir, const Vector3f& color, SG* out_SGs, uint64_t num_SGs, float sample_idx, float* lobe_weights, bool nonnegative);

} // namespace civet

#endif // CIVET_SG_H
