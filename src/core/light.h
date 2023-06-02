#ifndef CIVET_LIGHT_H
#define CIVET_LIGHT_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/medium.h>
#include <core/node.h>
#include <core/shader.h>
#include <core/interaction.h>

namespace civet {

enum class LightFlags : int {
	DeltaPosition = 1,
	DeltaDirection = 2,
	Area = 4,
	Infinite = 8
};

inline bool isDeltaLight(int flags) {
	return flags & (int)LightFlags::DeltaPosition ||
			flags & (int)LightFlags::DeltaDirection;
}

class Light {
public:
	Light(int flags, const Transform& ltw, const MediumInterface& mi, int n_samples = 1) :
			flags(flags), n_samples(std::max(1, n_samples)), light_to_world(ltw), world_to_light(inverse(ltw)) {
	}

	virtual Spectrum sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, float* pdf, VisibilityTester* vis) const = 0;
	virtual float pdf_Li(const Interaction& ref, const Vector3f& wi) const  = 0;
	virtual Spectrum power() const = 0;
	virtual void preprocess(const Scene& scene) {}

	const int flags;
	const int n_samples;
	const MediumInterface medium_interface;

protected:
	const Transform light_to_world, world_to_light;
};

class VisibilityTester {
public:
	VisibilityTester() {}
	VisibilityTester(const Interaction& p0, const Interaction& p1) :
			p0(p0), p1(p1) {}

	const Interaction &P0() const { return p0; }
	const Interaction &P1() const { return p1; }

	bool unoccluded(const Scene& scene) const;
	Spectrum Tr(const Scene& scene, Sampler& sampler) const;

private:
	Interaction p0, p1;
};

struct AttenuationFactor {
	float constant = 1.f;
	float linear = 0.14f;
	float quadratic = 0.07f;
};

class GLLight : public Node {
public:
	GLLight(const std::string& name, const NodeType type) :
			Node(name, type) {}

	virtual void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) = 0;

	virtual void init() = 0;

	unsigned int shadow_map;
	bool active;
	bool cast_shadow;

	Vector3f color;
	float power;

	bool should_update; ///< true only if need to regenerate shadow map, e.g. position or direction changes

protected:
	unsigned int FBO;
	unsigned int resolution;
};

class GLDirectionalLight : public GLLight {
public:
	GLDirectionalLight(const std::string& name, bool cascades = true, unsigned int res = 2048) :
			use_cascaded_shadows(cascades), GLLight(name, NodeType::DirectionalLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		power = 3.f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	GLDirectionalLight(const std::string& name, Vector3f dir, bool cascades = true, unsigned int res = 2048) :
			direction(dir), use_cascaded_shadows(cascades), GLLight(name, NodeType::DirectionalLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		power = 3.f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	~GLDirectionalLight();

	void generateShadowMap(Shader& shader);
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Vector3f direction{ 1.f, 0.f, 0.f };
	std::vector<Transform> light_space_mat{ 4 };

	unsigned int UBO;
	std::vector<float> cascade_levels;
	float frustum_fitting_factor = 100.f;
	bool use_cascaded_shadows;

	void init() override;

private:
	std::vector<Point3f> getFrustumCornersInWorldSpace(const Transform& projection, const Transform& view);
	Transform getLightSpaceMatrix(const float near_plane, const float far_plane);
	std::vector<Transform> getLightSpaceMatrices();
};

class GLPointLight : public GLLight {
public:
	GLPointLight(const std::string& name, unsigned int res = 1024) :
			GLLight(name, NodeType::PointLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		power = 100.f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	GLPointLight(const std::string& name, Point3f pos, unsigned int res = 1024) :
			position(pos), GLLight(name, NodeType::PointLight) {
		resolution = res;
		color = Vector3f(1, 1, 1);
		power = 100.f;
		active = true;
		cast_shadow = true;
		should_update = true;
	}

	~GLPointLight();

	void generateShadowMap(Shader& shader, float near_plane, float far_plane);
	void bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) override;

	Point3f position;
	float radius;
	AttenuationFactor attenuation;

	void init() override;
};

} // namespace civet

#endif // CIVET_LIGHT_H
