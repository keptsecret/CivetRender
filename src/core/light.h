#ifndef CIVET_LIGHT_H
#define CIVET_LIGHT_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/shader.h>

namespace civet {

class GLLight {
public:
	virtual void generateShadowMap(Shader& shader, float near_plane, float far_plane) = 0;

	virtual void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) = 0;

	unsigned int shadow_map;
	bool active;
	bool cast_shadow;

	Vector3f color;
	float intensity;

	virtual void init() = 0;

protected:
	unsigned int FBO;
	unsigned int resolution;

	bool should_update;
};

class GLDirectionalLight : public GLLight {
public:
	GLDirectionalLight(unsigned int res = 2048) {
		resolution = res;
		direction = Vector3f(1, 0, 0);
		color = Vector3f(1, 1, 1);
		intensity = 0.5f;
		active = true;
		cast_shadow = true;
	}

	GLDirectionalLight(Vector3f dir, unsigned int res = 2048) :
			direction(dir) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		intensity = 0.5f;
		active = true;
		cast_shadow = true;
	}

	void generateShadowMap(civet::Shader& shader, float near_plane, float far_plane) override;
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Vector3f direction;
	Transform light_space_mat;

	void init() override;
};

class GLPointLight : public GLLight {
public:
	GLPointLight(unsigned int res = 2048) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		intensity = 0.8f;
		active = true;
		cast_shadow = true;
	}

	GLPointLight(Point3f pos, unsigned int res = 2048) :
			position(pos) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		intensity = 0.8f;
		active = true;
		cast_shadow = true;
	}

	void generateShadowMap(civet::Shader& shader, float near_plane, float far_plane) override;
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Point3f position;

	void init() override;
};

} // namespace civet

#endif // CIVET_LIGHT_H
