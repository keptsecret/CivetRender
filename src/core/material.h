#ifndef CIVET_MATERIAL_H
#define CIVET_MATERIAL_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/shader.h>
#include <utils/memory.h>

namespace civet {

enum class TransportMode { Radiance,
	Importance };

class Material {
public:
	virtual ~Material();

	virtual void computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const = 0;

	static void bump(const std::shared_ptr<Texture<float>>& map, SurfaceInteraction* si);
};

struct GLTexture {
	unsigned int id;
	std::string type;
	std::string path;
};

class GLMaterial {
public:
	GLMaterial() {}

	void bind(Shader& shader, unsigned int tex_offset);
	void unbind(unsigned int tex_offset);

	Vector3f albedo{0.8f, 0.8f, 0.8f};
	float metallic = 0.f;
	float roughness = 0.5f;
	float ambient = 1.f;

	bool use_albedo_map = false;
	bool use_metallic_map = false;
	bool use_roughness_map = false;
	bool is_glossy_rough = false;
	bool use_ao_map = false;
	bool use_normal_map = false;
	bool use_bump_map = false;
	float bump_scale = 0.1f;
	std::vector<std::shared_ptr<GLTexture>> textures;
};

} // namespace civet

#endif // CIVET_MATERIAL_H
