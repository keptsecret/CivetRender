#ifndef CIVET_SAMPLER_H
#define CIVET_SAMPLER_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <utils/rng.h>

namespace civet {

// Utility functions ?
template <typename T>
void shuffle(T* sample, int count, int n_dims, RNG& rng);

class Sampler {
public:
	Sampler(int64_t spp);
	virtual ~Sampler() {}

	virtual void startPixel(const Point2i& p);

	virtual float get1D() = 0;
	virtual Point2f get2D() = 0;

	CameraSample getCameraSample(const Point2i& p_raster);
	void request1DArray(int n);
	void request2DArray(int n);
	virtual int roundCount(int n) const { return n; }

	const float* get1DArray(int n);
	const Point2f* get2DArray(int n);

	virtual bool startNextSample();

	virtual std::unique_ptr<Sampler> clone(int seed) = 0;
	virtual bool setSampleNumber(int64_t sample_num);

	const int64_t samples_per_pixel;

protected:
	Point2i current_pixel;
	int64_t current_pixel_sample_index;

	std::vector<int> samples1D_array_sizes, samples2D_array_sizes;
	std::vector<std::vector<float>> sample_array1D;
	std::vector<std::vector<Point2f>> sample_array2D;

private:
	size_t array1D_offset, array2D_offset;
};

class PixelSampler : public Sampler {
public:
	PixelSampler(int64_t spp, int n_sampled_dims);

	bool startNextSample();
	bool setSampleNumber(int64_t sample_num);

	float get1D();
	Point2f get2D();

protected:
	std::vector<std::vector<float>> samples1D;
	std::vector<std::vector<Point2f>> samples2D;
	int current1D_dimension = 0, current2D_dimension = 0;
	RNG rng;
};

class GlobalSampler : public Sampler {
public:
	GlobalSampler(int64_t spp);

	virtual int64_t getIndexForSample(int64_t sample_num) const = 0;
	virtual float sampleDimension(int64_t index, int dimension) const = 0;

	void startPixel(const Point2i &p) override;
	bool startNextSample() override;
	bool setSampleNumber(int64_t sample_num) override;

	float get1D() override;
	Point2f get2D() override;
private:
	int dimension;
	int64_t interval_sample_index;
	static const int array_start_dim = 5;
	int array_end_dim;
};

} // namespace civet

#endif // CIVET_SAMPLER_H
