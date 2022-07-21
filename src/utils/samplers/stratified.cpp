#include <utils/samplers/stratified.h>

namespace civet {

void stratifiedSample1D(float* sample, int n_samples, RNG& rng, bool jitter = true) {
	float inv_n_samples = 1.0f / n_samples;
	for (int i = 0; i < n_samples; i++) {
		float delta = jitter ? rng.uniformFloat() : 0.5f;
		sample[i] = std::min((i + delta) * inv_n_samples, OneMinusEpsilon);
	}
}

void stratifiedSample2D(Point2f* sample, int nx, int ny, RNG& rng, bool jitter = true) {
	float dx = 1.0f / nx, dy = 1.0f / ny;
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			float jx = jitter ? rng.uniformFloat() : 0.5f;
			float jy = jitter ? rng.uniformFloat() : 0.5f;
			sample->x = std::min((x + jx) * dx, OneMinusEpsilon);
			sample->y = std::min((y + jy) * dy, OneMinusEpsilon);
			sample++;
		}
	}
}

void latinHypercube(float* samples, int n_samples, int n_dim, RNG& rng) {
	float inv_n_samples = 1.0f / n_samples;
	for (int i = 0; i < n_samples; i++) {
		for (int j = 0; j < n_dim; j++) {
			float sj = (i + rng.uniformFloat()) * inv_n_samples;
			samples[n_dim * i + j] = std::min(sj, OneMinusEpsilon);
		}
	}

	for (int i = 0; i < n_dim; i++) {
		for (int j = 0; j < n_samples; j++) {
			int other = j + rng.uniformUInt32(n_samples - j);
			swapElem(samples[n_dim * j + i], samples[n_dim * other + i]);
		}
	}
}

void StratifiedSampler::startPixel(const Point2i& p) {
	for (size_t i = 0; i < samples1D.size(); i++) {
		stratifiedSample1D(&samples1D[i][0], x_pixel_samples * y_pixel_samples, rng, jitter_samples);
		shuffle(&samples1D[i][0], x_pixel_samples * y_pixel_samples, 1, rng);
	}
	for (size_t i = 0; i < samples2D.size(); i++) {
		stratifiedSample2D(&samples2D[i][0], x_pixel_samples, y_pixel_samples, rng, jitter_samples);
		shuffle(&samples2D[i][0], x_pixel_samples * y_pixel_samples, 1, rng);
	}

	for (size_t i = 0; i < samples1D_array_sizes.size(); i++) {
		for (int64_t j = 0; j < samples_per_pixel; j++) {
			int count = samples1D_array_sizes[i];
			stratifiedSample1D(&samples1D[i][j * count], count, rng, jitter_samples);
			shuffle(&samples1D[i][j * count], count, 1, rng);
		}
	}
	for (size_t i = 0; i < samples2D_array_sizes.size(); i++) {
		for (int64_t j = 0; j < samples_per_pixel; j++) {
			int count = samples2D_array_sizes[i];
			latinHypercube(&samples2D[i][j * count].x, count, 2, rng);
		}
	}

	PixelSampler::startPixel(p);
}

std::unique_ptr<Sampler> StratifiedSampler::clone(int seed) {
	StratifiedSampler *ss = new StratifiedSampler(*this);
	ss->rng.setSequence(seed);
	return std::unique_ptr<Sampler>(ss);
}

} // namespace civet