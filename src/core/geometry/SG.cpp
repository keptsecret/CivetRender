#include <core/geometry/SG.h>

#include <utils/sampling.h>

namespace civet {

static const int MaxSGCount = 12;
static SG defaultInitialGuess[MaxSGCount];

void generateUniformSGs(SG* out_SGs, uint64_t num_SGs, SGDistribution distribution) {
	const uint64_t N = distribution == SGDistribution::Hemispherical ? num_SGs : 2 * num_SGs;
	Vector3f means[MaxSGCount * 2];

	float inc = Pi * (3.f - std::sqrt(5.f));
	float off = 2.0f / N;
	for (uint64_t k = 0; k < N; k++) {
		float y = k * off - 1.0f + (off / 2.0f);
		float r = std::sqrt(1.0f - y * y);
		float phi = k * inc;
		means[k] = Vector3f(std::cos(phi) * r, std::sin(phi) * r, y);
	}

	uint64_t currSG = 0;
	for (uint64_t i = 0; i < N; i++) {
		// for the hemisphere only accept points on correct side
		if (distribution == SGDistribution::Spherical || dot(means[i], Vector3f(0.0f, 0.0f, 1.0f)) >= 0.0f) {
			SG sample;
			sample.axis = normalize(means[i]);
			out_SGs[currSG++] = sample;
		}
	}

	float min_DP = 1.0f;
	for (uint64_t i = 1; i < num_SGs; ++i) {
		Vector3f h = normalize(out_SGs[i].axis + out_SGs[0].axis);
		min_DP = std::min(min_DP, dot(h, out_SGs[0].axis));
	}

	float sharpness = (std::log(0.65f) * num_SGs) / (min_DP - 1.0f);

	for (uint32_t i = 0; i < num_SGs; ++i) {
		out_SGs[i].sharpness = sharpness;
	}

	const uint64_t sampleCount = 2048;
	Vector2f samples[sampleCount];
	generateHammersleySamples2D(samples, sampleCount);

	for (uint32_t i = 0; i < num_SGs; ++i) {
		out_SGs[i].basis_sq_integral_over_domain = 0.0f;
	}

	for (uint64_t i = 0; i < sampleCount; ++i) {
		Point2f sample(samples[i]);
		Vector3f dir = distribution == SGDistribution::Hemispherical ? uniformSampleHemisphere(sample) : uniformSampleSphere(sample);
		for (uint32_t j = 0; j < num_SGs; ++j) {
			float weight = std::exp(out_SGs[j].sharpness * (dot(dir, out_SGs[j].axis) - 1.0f));
			out_SGs[j].basis_sq_integral_over_domain += (weight * weight - out_SGs[j].basis_sq_integral_over_domain) / float(i + 1);
		}
	}
}

void initializeSGSolver(uint64_t num_SGs, SGDistribution distribution) {
	generateUniformSGs(defaultInitialGuess, num_SGs, distribution);
}

const SG* initialGuess() {
	return defaultInitialGuess;
}

void projectOntoSGs(const Vector3f& dir, const Vector3f& color, SG* out_SGs, uint64_t num_SGs) {
	for (uint64_t i = 0; i < num_SGs; i++) {
		SG sg1, sg2;
		sg1.amplitude = out_SGs[i].amplitude;
		sg1.axis = out_SGs[i].axis;
		sg1.sharpness = out_SGs[i].sharpness;
		sg2.amplitude = color;
		sg2.axis = normalize(dir);

		if (dot(dir, sg1.axis) > 0.f) {
			float dot12 = dot(sg1.axis, sg2.axis);
			float factor = (dot12 - 1.f) * sg1.sharpness;
			float wgt = std::exp(factor);
			out_SGs[i].amplitude += sg2.amplitude * wgt;
		}
	}
}

// Do a projection of the colors onto the SG's
static void solveProjection(SGSolveParam& params) {
	// Project color samples onto the SGs
	for (int i = 0; i < params.num_samples; i++) {
		projectOntoSGs(params.x_samples[i], params.y_samples[i], params.out_SGs, params.num_SGs);
	}

	// Weight the samples by the monte carlo factor for uniformly sampling the sphere
	float MC_factor = ((4.f * Pi) / params.num_samples);
	for (int i = 0; i < params.num_SGs; i++) {
		params.out_SGs[i].amplitude *= MC_factor;
	}
}

// Accumulates a single sample for computing a set of SG's using a running average. This technique and the code it's based
// on was provided by Thomas Roughton in the following article: http://torust.me/rendering/irradiance-caching/spherical-gaussians/2018/09/21/spherical-gaussians.html
void SGRunningAverage(const Vector3f& dir, const Vector3f& color, SG* out_SGs, uint64_t num_SGs, float sample_idx, float* lobe_weights, bool nonnegative) {
	float sample_weight_scale = 1.0f / (sample_idx + 1);

	float sample_lobe_weights[MaxSGCount] = {};
	Vector3f curr_estimate;

	for (uint64_t lobe_idx = 0; lobe_idx < num_SGs; lobe_idx++) {
		float dot_prod = dot(out_SGs[lobe_idx].axis, dir);
		float weight = exp(out_SGs[lobe_idx].sharpness * (dot_prod - 1.0f));
		curr_estimate += out_SGs[lobe_idx].amplitude * weight;

		sample_lobe_weights[lobe_idx] = weight;
	}

	for (uint64_t lobe_idx = 0; lobe_idx < num_SGs; lobe_idx++) {
		float weight = sample_lobe_weights[lobe_idx];
		if (weight == 0.0f) {
			continue;
		}

		float spherical_integral_guess = weight * weight;

		lobe_weights[lobe_idx] += (spherical_integral_guess - lobe_weights[lobe_idx]) * sample_weight_scale;

		// Clamp the spherical integral estimate to at least the true value to reduce variance.
		float spherical_integral = std::max(lobe_weights[lobe_idx], out_SGs[lobe_idx].basis_sq_integral_over_domain);

		Vector3f other_contributions = curr_estimate - out_SGs[lobe_idx].amplitude * weight;
		Vector3f newValue = (color - other_contributions) * (weight / spherical_integral);

		out_SGs[lobe_idx].amplitude += (newValue - out_SGs[lobe_idx].amplitude) * sample_weight_scale;

		if (nonnegative) {
			out_SGs[lobe_idx].amplitude.x = std::max(out_SGs[lobe_idx].amplitude.x, 0.0f);
			out_SGs[lobe_idx].amplitude.y = std::max(out_SGs[lobe_idx].amplitude.y, 0.0f);
			out_SGs[lobe_idx].amplitude.z = std::max(out_SGs[lobe_idx].amplitude.z, 0.0f);
		}
	}
}

static void solveRunningAverage(SGSolveParam& params, bool nonnegative) {
	float lobeWeights[MaxSGCount] = {};

	// Project color samples onto the SGs
	for (uint32_t i = 0; i < params.num_samples; i++) {
		SGRunningAverage(params.x_samples[i], params.y_samples[i], params.out_SGs, params.num_SGs, (float)i, lobeWeights, nonnegative);
	}
}

void solveSGs(SGSolveParam& params) {
	for (uint64_t i = 0; i < params.num_SGs; i++) {
		params.out_SGs[i] = defaultInitialGuess[i];
	}

	solveRunningAverage(params, true);
}


} // namespace civet