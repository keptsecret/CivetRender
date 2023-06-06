#ifndef CIVET_ILLUMINANCEFIELD_H
#define CIVET_ILLUMINANCEFIELD_H

#include <core/shader.h>
#include <core/civet.h>
#include <core/geometry/SG.h>
#include <core/geometry/vecmath.h>
#include <utils/sampling.h>

namespace civet {

class Editor;

struct SGProbe {
	static const uint64_t SG_count = 12; // each probe is represented by 12 SGs
	uint64_t num_samples = 0;
	uint64_t curr_sample = 0;
	std::vector<Vector3f> sample_dirs;
	std::vector<Vector3f> samples;

	SG projected_result[SG_count];
	float running_avg_weights[SG_count] = {};

	std::vector<Vector3f> amplitudes; // each SG has different amplitudes determined by direction

	void init(uint64_t ns) {
		curr_sample = 0;
		num_samples = ns;

		if (num_samples != sample_dirs.size()) {
			sample_dirs.clear();
			sample_dirs.resize(num_samples);
		}
		if (num_samples != samples.size()) {
			samples.clear();
			samples.resize(num_samples);
		}

		for (uint64_t i = 0; i < SG_count; i++) {
			projected_result[i].amplitude = Vector3f(); ///< possibly redundant
			running_avg_weights[i] = 0.f;
		}

		amplitudes.resize(SG_count);
	}

	Vector3f sampleDirection(const Point2f& sample) { // TODO: check this, assuming we sample for SG on all directions instead of hemisphere
		return uniformSampleSphere(sample);
	}

	void addSample(const Vector3f& sample_dir, uint64_t sample_idx, const Vector3f& sample) {
		sample_dirs[curr_sample] = sample_dir;
		samples[curr_sample] = sample;
		curr_sample++;

//		SGRunningAverage(sample_dir, sample, projected_result, SG_count, (float)sample_idx, running_avg_weights, true);
		projectOntoSGs(sample_dir, sample, projected_result, SG_count);
	}

	void bakeResult() {
		SG lobes[SG_count];

		SGSolveParam params;
		params.num_SGs = SG_count;
		params.out_SGs = lobes;
		params.x_samples = sample_dirs.data();
		params.y_samples = samples.data();
		params.num_samples = num_samples;
		solveSGs(params);

		for (int i = 0; i < SG_count; i++) {
			amplitudes[i].x = clamp(lobes[i].amplitude.x, 0, Infinity);
			amplitudes[i].y = clamp(lobes[i].amplitude.y, 0, Infinity);
			amplitudes[i].z = clamp(lobes[i].amplitude.z, 0, Infinity);
		}
	}
};

class IlluminanceField {
public:
	~IlluminanceField();

	void initialize();
	void bake(const Scene& scene);
	void bind(Shader& shader, int tex_offset);

	void fitGridToBounds(const Bounds3f& bounds);

	bool hasBakeData() const { return has_bake_data; }

private:
	friend Editor;
	Point3i probe_grid_size{ 16, 16, 16 };
	float cell_size = 10.f;
	uint64_t rays_per_probe = 1024;
	Point3f corner_position{ 0, 0, 0 };

	std::vector<SGProbe> probes;
	Vector3f SG_directions[SGProbe::SG_count]; // directions are same for each probe
	float SG_sharpness; // sharpness is same for all initialized SGs

	unsigned int SG_data_texture_array;
	bool has_bake_data = false;
};

} // namespace civet

#endif // CIVET_ILLUMINANCEFIELD_H
