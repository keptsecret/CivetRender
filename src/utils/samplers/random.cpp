#include <utils/samplers/random.h>

namespace civet {

void RandomSampler::startPixel(const Point2i& p) {
	for (size_t i = 0; i < sample_array1D.size(); i++) {
		for (size_t j = 0; j < sample_array1D[i].size(); j++) {
			sample_array1D[i][j] = rng.uniformFloat();
		}
	}

	for (size_t i = 0; i < sample_array2D.size(); i++) {
		for (size_t j = 0; j < sample_array2D[i].size(); j++) {
			sample_array2D[i][j] = { rng.uniformFloat(), rng.uniformFloat() };
		}
	}
	Sampler::startPixel(p);
}

float RandomSampler::get1D() {
	return rng.uniformFloat();
}

Point2f RandomSampler::get2D() {
	return { rng.uniformFloat(), rng.uniformFloat() };
}

std::unique_ptr<Sampler> RandomSampler::clone(int seed) {
	RandomSampler* rs = new RandomSampler(*this);
	rs->rng.setSequence(seed);
	return std::unique_ptr<Sampler>(rs);
}

} // namespace civet