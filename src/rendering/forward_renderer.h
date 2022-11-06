#ifndef CIVET_FORWARD_RENDERER_H
#define CIVET_FORWARD_RENDERER_H

#include <core/camera.h>
#include <core/civet.h>
#include <core/light.h>
#include <core/mesh.h>
#include <shaders/simple_forward.h>

namespace civet {

class ForwardRenderer {
public:
	ForwardRenderer(unsigned int w, unsigned int h) :
			width(w), height(h) {}

	void draw(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights,
			SimpleForwardShader* forward, const unsigned int shadow_res, Shader* dir_shadow = nullptr, Shader* omni_shadow = nullptr);

	void setModelMat(const Transform& model) { model_mat = model; }
	void setViewMat(const Transform& view) { view_mat = view; }
	void setProjectionMat(const Transform& projection) { projection_mat = projection; }
	void setCamera(GLCamera* cam) { camera = cam; }

private:
	void generateShadowMaps(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights,
			const unsigned int shadow_res, Shader* depth_shader = nullptr, Shader* depth_cube_shader = nullptr);
	void renderScene(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights, SimpleForwardShader* shader);

	unsigned int width, height;
	Transform model_mat, view_mat, projection_mat;
	GLCamera* camera;
};

} // namespace civet

#endif // CIVET_FORWARD_RENDERER_H
