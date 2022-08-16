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

private:
	virtual void init() = 0;

protected:
	unsigned int FBO;
	unsigned int resolution;

	bool active;
	bool cast_shadow;
};

class GLDirectionalLight : public GLLight {
public:
	GLDirectionalLight(unsigned int res) {
		resolution = res;
	}

	GLDirectionalLight(unsigned int res, Vector3f dir) :
			direction(dir) {
		resolution = res;
	}

	void generateShadowMap(civet::Shader& shader, float near_plane, float far_plane) override;
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Vector3f direction;

private:
	void init() override;
};

class GLPointLight : public GLLight {
public:
	//TODO: needs fixing - cube shadow maps currently output black in testing
	GLPointLight(unsigned int res) {
		resolution = res;
	}

	GLPointLight(unsigned int res, Point3f pos) :
			position(pos) {
		resolution = res;
	}

	void generateShadowMap(civet::Shader& shader, float near_plane, float far_plane) override;
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Point3f position;

private:
	void init() override;
};

} // namespace civet

#endif // CIVET_LIGHT_H
