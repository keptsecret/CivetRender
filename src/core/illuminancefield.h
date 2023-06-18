#ifndef CIVET_ILLUMINANCEFIELD_H
#define CIVET_ILLUMINANCEFIELD_H

#include <core/civet.h>
#include <core/geometry/SG.h>
#include <core/geometry/vecmath.h>
#include <core/shader.h>
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

		SGRunningAverage(sample_dir, sample, projected_result, SG_count, (float)sample_idx, running_avg_weights, true);
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

	void testPathtracer(const Scene& scene);

private:
	friend Editor;
	Point3f corner_position{ 0, 0, 0 };
	Point3i probe_grid_size{ 16, 16, 16 };
	Vector3f cell_dim{ 1, 1, 1 };
	uint64_t rays_per_probe = 4096;

	std::vector<SGProbe> probes;
	Vector3f SG_directions[SGProbe::SG_count];	// directions are same for each probe
	float SG_sharpness;	// sharpness is same for all initialized SGs

	uint64_t num_distance_samples = 128;
	float distance_lobe_size = 0.5f;

	unsigned int FBO;
	int cubemap_resolution = 256;
	int octahedral_resolution = 128;
	unsigned int distance_cubemap;

	unsigned int SG_data_texture_array;
	unsigned int filtered_distance_texture_array;
	unsigned int sphere_samples_UBO;

	Shader distance_cubemap_shader;
	Shader background_distance_shader;
	Shader filtered_distance_shader;
	bool has_bake_data = false;
};

} // namespace civet

#endif // CIVET_ILLUMINANCEFIELD_H
