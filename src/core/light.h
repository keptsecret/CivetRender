#ifndef CIVET_LIGHT_H
#define CIVET_LIGHT_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/shader.h>
#include <core/node.h>

namespace civet {

struct AttenuationFactor {
	float constant = 1.f;
	float linear = 0.14f;
	float quadratic = 0.07f;
};

class GLLight : public Node {
public:
	GLLight(const std::string& name, const NodeType type) :
			Node(name, type) {}

	virtual void generateShadowMap(Shader& shader, float near_plane, float far_plane) = 0;
	virtual void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) = 0;

	virtual void init() = 0;

	unsigned int shadow_map;
	bool active;
	bool cast_shadow;

	Vector3f color;
	float intensity;

	bool should_update;	///< true only if need to regenerate shadow map, e.g. position or direction changes

protected:
	unsigned int FBO;
	unsigned int resolution;
};

class GLDirectionalLight : public GLLight {
public:
	GLDirectionalLight(const std::string& name, unsigned int res = 2048) :
			GLLight(name, DirectionalLight) {
		resolution = res;
		direction = Vector3f(1, 0, 0);
		color = Vector3f(1, 1, 1);
		intensity = 0.5f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	GLDirectionalLight(const std::string& name, Vector3f dir, unsigned int res = 2048) :
			direction(dir), GLLight(name, DirectionalLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		intensity = 0.5f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	void generateShadowMap(civet::Shader& shader, float near_plane, float far_plane) override;
	void generateShadowMap(civet::Shader& shader, Bounds3f frustum);
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Vector3f direction;
	Transform light_space_mat;

	void init() override;
};

class GLPointLight : public GLLight {
public:
	GLPointLight(const std::string& name, unsigned int res = 2048) :
			GLLight(name, PointLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		intensity = 0.8f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	GLPointLight(const std::string& name, Point3f pos, unsigned int res = 2048) :
			position(pos), GLLight(name, PointLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		intensity = 0.8f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	void generateShadowMap(civet::Shader& shader, float near_plane, float far_plane) override;
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Point3f position;
	float far_plane;
	AttenuationFactor attenuation;

	void init() override;
};

} // namespace civet

#endif // CIVET_LIGHT_H
