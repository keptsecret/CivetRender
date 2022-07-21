#include <core/sampler.h>

#include <core/camera.h>

namespace civet {

template <typename T>
void shuffle(T* sample, int count, int n_dims, RNG& rng) {
	for (int i = 0; i < count; i++) {
		int other = i + rng.uniformUInt32(count - i);
		for (int j = 0; j < n_dims; j++) {
			swapElem(sample[n_dims * i + j], sample[n_dims * other + j]);
		}
	}
}

Sampler::Sampler(int64_t spp) :
		samples_per_pixel(spp) {
}

void Sampler::startPixel(const Point2i& p) {
	current_pixel = p;
	current_pixel_sample_index = 0;
	array1D_offset = array2D_offset = 0;
}

CameraSample Sampler::getCameraSample(const Point2i& p_raster) {
	CameraSample cs;
	cs.p_film = Point2f(p_raster.x, p_raster.y) + get2D();
	cs.time = get1D();
	cs.p_lens = get2D();
	return cs;
}

void Sampler::request1DArray(int n) {
	samples1D_array_sizes.push_back(n);
	sample_array1D.push_back(std::vector<float>(n * samples_per_pixel));
}

void Sampler::request2DArray(int n) {
	samples2D_array_sizes.push_back(n);
	sample_array2D.push_back(std::vector<Point2f>(n * samples_per_pixel));
}

const float* Sampler::get1DArray(int n) {
	if (array1D_offset == sample_array1D.size()) {
		return nullptr;
	}
	return &sample_array1D[array1D_offset++][current_pixel_sample_index * n];
}

const Point2f* Sampler::get2DArray(int n) {
	if (array2D_offset == sample_array2D.size()) {
		return nullptr;
	}
	return &sample_array2D[array2D_offset++][current_pixel_sample_index * n];
}

bool Sampler::startNextSample() {
	array1D_offset = array2D_offset = 0;
	return ++current_pixel_sample_index < samples_per_pixel;
}

bool Sampler::setSampleNumber(int64_t sample_num) {
	array1D_offset = array2D_offset = 0;
	current_pixel_sample_index = sample_num;
	return current_pixel_sample_index < samples_per_pixel;
}

PixelSampler::PixelSampler(int64_t spp, int n_sampled_dims) :
		Sampler(spp) {
	for (int i = 0; i < n_sampled_dims; i++) {
		samples1D.push_back(std::vector<float>(spp));
		samples2D.push_back(std::vector<Point2f>(spp));
	}
}

bool PixelSampler::startNextSample() {
	current1D_dimension = current2D_dimension = 0;
	return Sampler::startNextSample();
}

bool PixelSampler::setSampleNumber(int64_t sample_num) {
	current1D_dimension = current2D_dimension = 0;
	return Sampler::setSampleNumber(sample_num);
}

float PixelSampler::get1D() {
	if (current1D_dimension < samples1D.size()) {
		return samples1D[current1D_dimension++][current_pixel_sample_index];
	} else {
		return rng.uniformFloat();
	}
}

Point2f PixelSampler::get2D() {
	if (current2D_dimension < samples2D.size()) {
		return samples2D[current2D_dimension++][current_pixel_sample_index];
	} else {
		return Point2f(rng.uniformFloat(), rng.uniformFloat());
	}
}

GlobalSampler::GlobalSampler(int64_t spp) :
		Sampler(spp) {
}

void GlobalSampler::startPixel(const Point2i& p) {
	Sampler::startPixel(p);
	dimension = 0;
	interval_sample_index = getIndexForSample(0);

	array_end_dim = array_start_dim + sample_array1D.size() + 2 * sample_array2D.size();

	for (size_t i = 0; i < samples1D_array_sizes.size(); i++) {
		int n_samples = samples1D_array_sizes[i] * samples_per_pixel;
		for (int j = 0; j < n_samples; j++) {
			int64_t index = getIndexForSample(j);
			sample_array1D[i][j] = sampleDimension(index, array_start_dim + i);
		}
	}

	int dim = array_start_dim + samples1D_array_sizes.size();
	for (size_t i = 0; i < samples2D_array_sizes.size(); i++) {
		int n_samples = samples2D_array_sizes[i] * samples_per_pixel;
		for (int j = 0; j < n_samples; j++) {
			int64_t index = getIndexForSample(j);
			sample_array2D[i][j].x = sampleDimension(index, dim);
			sample_array2D[i][j].y = sampleDimension(index, dim + 1);
		}
		dim += 2;
	}
}

bool GlobalSampler::startNextSample() {
	dimension = 0;
	interval_sample_index = getIndexForSample(current_pixel_sample_index + 1);
	return Sampler::startNextSample();
}

bool GlobalSampler::setSampleNumber(int64_t sample_num) {
	dimension = 0;
	interval_sample_index = getIndexForSample(sample_num);
	return Sampler::setSampleNumber(sample_num);
}

float GlobalSampler::get1D() {
	if (dimension >= array_start_dim && dimension < array_end_dim) {
		dimension = array_end_dim;
	}
	return sampleDimension(interval_sample_index, dimension++);
}

Point2f GlobalSampler::get2D() {
	if (dimension + 1 >= array_start_dim && dimension < array_end_dim) {
		dimension = array_end_dim;
	}
	Point2f p(sampleDimension(interval_sample_index, dimension), sampleDimension(interval_sample_index, dimension + 1));
	dimension += 2;
	return p;
}

} // namespace civet