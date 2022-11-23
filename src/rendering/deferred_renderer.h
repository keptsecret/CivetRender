#ifndef CIVET_DEFERRED_RENDERER_H
#define CIVET_DEFERRED_RENDERER_H

#include <core/camera.h>
#include <core/civet.h>
#include <core/light.h>
#include <core/mesh.h>
#include <rendering/gbuffer.h>

namespace civet {

class DeferredRenderer {
protected:
	DeferredRenderer() {}

public:
	static DeferredRenderer* getSingleton();

	void init(unsigned int w, unsigned int h);

	void draw(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights);

	void setModelMat(const Transform& model) { model_mat = model; }
	void setViewMat(const Transform& view) { view_mat = view; }
	void setProjectionMat(const Transform& projection) { projection_mat = projection; }
	void setCamera(GLCamera* cam) { camera = cam; }

private:
	void geometryPass(GLModel& model);
	void lightsPass(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights);

	void pointLightsPass(GLModel& model, std::vector<GLPointLight>& lights);
	void dirLightsPass(GLModel& model, std::vector<GLDirectionalLight>& lights);

	float getBoundingSphere(GLPointLight& light);

	unsigned int width, height;
	Transform model_mat, view_mat, projection_mat;
	GLCamera* camera;
	GBuffer gbuffer;

	Shader geometry_pass_shader;
	Shader pointlight_pass_shader;
	Shader dirlight_pass_shader;

	GLModel bounding_sphere;
	GLModel bounding_quad;
};

} // namespace civet

#endif // CIVET_DEFERRED_RENDERER_H
