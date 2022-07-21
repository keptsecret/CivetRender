#ifndef CIVET_STRATIFIED_H
#define CIVET_STRATIFIED_H

#include <core/sampler.h>

namespace civet {

class StratifiedSampler : public PixelSampler {
public:
	StratifiedSampler(int x_samples, int y_samples, bool jitter, int n_sampled_dims) :
			PixelSampler(x_samples * y_samples, n_sampled_dims), x_pixel_samples(x_samples), y_pixel_samples(y_samples), jitter_samples(jitter) {}

	void startPixel(const Point2i& p);

	std::unique_ptr<Sampler> clone(int seed);

private:
	const int x_pixel_samples, y_pixel_samples;
	const bool jitter_samples;
};

} // namespace civet

#endif // CIVET_STRATIFIED_H
