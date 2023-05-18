#include <rendering/forward_renderer.h>

namespace civet {
void ForwardRenderer::draw(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights,
		SimpleForwardShader* forward, const unsigned int shadow_res, Shader* dir_shadow, Shader* omni_shadow) {
	glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	generateShadowMaps(model, dir_lights, point_lights, shadow_res, dir_shadow, omni_shadow);
	renderScene(model, dir_lights, point_lights, forward);
}

void ForwardRenderer::generateShadowMaps(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights,
		const unsigned int shadow_res, Shader* depth_shader, Shader* depth_cube_shader) {
	if (depth_shader) {
		// Render depth map
		float near_plane = -5.0f, far_plane = 5.0f;

		depth_shader->use();
		// loop through lights and generate shadow maps
		glViewport(0, 0, shadow_res, shadow_res);
		glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
		for (auto& light : dir_lights) {
			if (light.should_update) {
				depth_shader->setMat4("model", model_mat.m);
				light.generateShadowMap(*depth_shader);
				model.draw(*depth_shader, 2); ///< change tex_offset possibly
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
			}
		}
	}

	if (depth_cube_shader) {
		// Render depth cube map
		float near = 1.0f, far = 25.0f;

		depth_cube_shader->use();
		// loop through lights and generate shadow maps for point lights
		glViewport(0, 0, shadow_res, shadow_res);
		glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
		for (auto& light : point_lights) {
			if (light.should_update) {
				depth_cube_shader->setMat4("model", model_mat.m);
				light.generateShadowMap(*depth_cube_shader, near, far);
				model.draw(*depth_cube_shader, 2);
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
			}
		}
	}
}

void ForwardRenderer::renderScene(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights,
		SimpleForwardShader* shader) {
	// Render the scene
	glViewport(0, 0, width, height);
	glCullFace(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	shader->use();
	shader->setMat4("projection", projection_mat.m);
	shader->setMat4("view", view_mat.m);
	shader->setMat4("model", model_mat.m);

	shader->setVec3("viewPos", Vector3f(camera->position));
	shader->setFloat("material.shininess", 64.0f);

	int tex_offset = 0;
	///< put this into bind lights function in shader, possibly
	{
		for (int i = 0; i < dir_lights.size(); i++) {
			shader->setMat4(("dirLights[" + std::to_string(i) + "].light_space_mat").c_str(), dir_lights[i].light_space_mat[0].m);
			dir_lights[i].bindShadowMap(*shader, ("dirLights[" + std::to_string(i) + "].shadow_map").c_str(), tex_offset++);
		}

		for (int i = 0; i < point_lights.size(); i++) {
			shader->setFloat(("pointLights[" + std::to_string(i) + "].far_plane").c_str(), point_lights[i].radius);
			point_lights[i].bindShadowMap(*shader, ("pointLights[" + std::to_string(i) + "].shadow_map").c_str(), tex_offset++);
		}
	}

	model.draw(*shader, tex_offset);
}

} // namespace civet