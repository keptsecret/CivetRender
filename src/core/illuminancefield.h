#ifndef CIVET_ILLUMINANCEFIELD_H
#define CIVET_ILLUMINANCEFIELD_H

#include <core/civet.h>
#include <core/geometry/SG.h>
#include <core/geometry/vecmath.h>
#include <core/shader.h>
#include <utils/sampling.h>

namespace civet {

class Editor;

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
	Vector3f mapToDirection(float x, float y, int s);

private:
	friend Editor;
	Point3f corner_position{ 0, 0, 0 };
	Point3i probe_grid_size{ 16, 16, 16 };
	Vector3f cell_dim{ 30, 40, 30 };
	uint64_t rays_per_texel_radiance = 1;

	uint64_t num_irradiance_samples = 4096;
	uint64_t num_distance_samples = 256;
	float irradiance_lobe_size = 0.99f;
	float distance_lobe_size = 0.5f;

	std::vector<std::vector<Vector3f>> probe_radiance;
	std::vector<std::vector<Vector2f>> probe_distance;

	unsigned int FBO;
	int cubemap_resolution = 256;
	int octahedral_resolution = 128;
	unsigned int radiance_cubemap;
	unsigned int distance_cubemap;
	std::vector<std::vector<Vector3f>> radiance_cubemap_data;
	std::vector<std::vector<Vector2f>> distance_cubemap_data;

	unsigned int radiance_texture_array;
	unsigned int distance_texture_array;
	unsigned int irradiance_texture_array;
	unsigned int filtered_distance_texture_array;
	unsigned int sphere_samples_UBO;

	Shader irradiance_shader;
	bool has_bake_data = false;
};

} // namespace civet

#endif // CIVET_ILLUMINANCEFIELD_H
