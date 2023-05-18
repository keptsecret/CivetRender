#ifndef CIVET_DEFERRED_RENDERER_H
#define CIVET_DEFERRED_RENDERER_H

#include <core/camera.h>
#include <core/civet.h>
#include <core/light.h>
#include <core/mesh.h>
#include <core/scene.h>
#include <rendering/gbuffer.h>

namespace civet {

class DeferredRenderer {
protected:
	DeferredRenderer() {}

public:
	static DeferredRenderer* getSingleton();

	void init(unsigned int w, unsigned int h);

	void draw(Scene& scene);

	void setViewMat(const Transform& view) { view_mat = view; }
	void setProjectionMat(const Transform& projection) { projection_mat = projection; }
	void setCamera(GLCamera* cam) { camera = cam; }

private:
	void geometryPass(GLModel& model);
	void lightsPass(GLModel& model, std::vector<std::shared_ptr<GLDirectionalLight>>& dir_lights, std::vector<std::shared_ptr<GLPointLight>>& point_lights);
	void postProcessPass(Scene& scene);
	void finalPass();

	void pointLightPass(GLModel& model, GLPointLight& light);
	void stencilPass(GLModel& model, GLPointLight& light);
	void dirLightPass(GLModel& model, GLDirectionalLight& light);

	float getBoundingSphere(GLPointLight& light);
	void generateShadowMaps(GLModel& model, std::vector<std::shared_ptr<GLDirectionalLight>>& dir_lights, std::vector<std::shared_ptr<GLPointLight>>& point_lights);

	unsigned int width, height;
	const unsigned int shadow_res = 4096;
	Transform view_mat, projection_mat;
	GLCamera* camera;
	GBuffer gbuffer;

	Shader geometry_pass_shader;
	Shader pointlight_pass_shader;
	Shader dirlight_pass_shader;
	Shader stencil_pass_shader; ///< doesn't actually do anything except population depth and stencil
	Shader postprocess_shader;

	Shader depth_shader;
	Shader depth_cascade_shader;
	Shader depth_cube_shader;

	GLModel bounding_sphere;
	GLModel bounding_quad;
};

} // namespace civet

#endif // CIVET_DEFERRED_RENDERER_H
