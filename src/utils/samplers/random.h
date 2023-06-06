#ifndef CIVET_RANDOM_H
#define CIVET_RANDOM_H

#include <core/sampler.h>
#include <utils/rng.h>

namespace civet {

class RandomSampler : public Sampler {
public:
	RandomSampler(int ns, int seed = 0) :
			Sampler(ns), rng(seed) {}
	void startPixel(const civet::Point2i& p) override;
	float get1D() override;
	Point2f get2D() override;
	std::unique_ptr<Sampler> clone(int seed) override;

private:
	RNG rng;
};

} // namespace civet

#endif // CIVET_RANDOM_H
