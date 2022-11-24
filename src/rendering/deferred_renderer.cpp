#include <rendering/deferred_renderer.h>

namespace civet {

void DeferredRenderer::init(unsigned int w, unsigned int h) {
	width = w, height = h;
	if (!gbuffer.init(width, height)) {
		printf("ERROR::DeferredRenderer: GBuffer init error");
	}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCheckError("ERROR::Engine::init: OpenGL error code");

	bounding_sphere = GLModel("../civet/resources/basic-meshes/sphere.obj");
	bounding_quad = GLModel("../civet/resources/basic-meshes/quad.obj");

	geometry_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_geometry_frag.glsl");
	pointlight_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_pointlight_pass_frag.glsl");
	dirlight_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_dirlight_pass_frag.glsl");
	stencil_pass_shader = Shader("../civet/src/shaders/null_vert.glsl", "../civet/src/shaders/null_frag.glsl");

	pointlight_pass_shader.use();
	pointlight_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	pointlight_pass_shader.setInt("DiffuseMap", GBuffer::GBUFFER_TEXTURE_DIFFUSE);
	pointlight_pass_shader.setInt("SpecularMap", GBuffer::GBUFFER_TEXTURE_SPECULAR);
	pointlight_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	dirlight_pass_shader.use();
	dirlight_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	dirlight_pass_shader.setInt("DiffuseMap", GBuffer::GBUFFER_TEXTURE_DIFFUSE);
	dirlight_pass_shader.setInt("SpecularMap", GBuffer::GBUFFER_TEXTURE_SPECULAR);
	dirlight_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	depth_shader = Shader("../civet/src/shaders/light_depth_vert.glsl", "../civet/src/shaders/light_depth_frag.glsl");
	depth_cube_shader = Shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/light_cube_depth_frag.glsl", "../civet/src/shaders/light_cube_depth_geom.glsl");
}

void DeferredRenderer::draw(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights) {
	gbuffer.start();
	generateShadowMaps(model, dir_lights, point_lights);

	glViewport(0, 0, width, height);
	glCullFace(GL_BACK);
	geometryPass(model);
	lightsPass(model, dir_lights, point_lights);
	finalPass();
}

void DeferredRenderer::geometryPass(GLModel& model) {
	geometry_pass_shader.use();
	gbuffer.bindGeomPass();

	glDepthMask(GL_TRUE); // only update depth buffer on geometry pass
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	geometry_pass_shader.use();
	geometry_pass_shader.setMat4("projection", projection_mat.m);
	geometry_pass_shader.setMat4("view", view_mat.m);
	geometry_pass_shader.setMat4("model", model_mat.m);
	model.draw(geometry_pass_shader, 0);

	glDepthMask(GL_FALSE);
}

void DeferredRenderer::lightsPass(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights) {
	glEnable(GL_STENCIL_TEST);	///< render light only if it passes test
	// handle point lights
	stencil_pass_shader.use();
	stencil_pass_shader.setMat4("projection", projection_mat.m);
	stencil_pass_shader.setMat4("view", view_mat.m);

	pointlight_pass_shader.use();
	pointlight_pass_shader.setMat4("projection", projection_mat.m);
	pointlight_pass_shader.setMat4("view", view_mat.m);

	pointlight_pass_shader.setVec3("viewPos", Vector3f(camera->position));
	pointlight_pass_shader.setVec2("screenSize", width, height);
	pointlight_pass_shader.setFloat("material.shininess", 64.0f);
	for (auto& light : point_lights) {
		stencilPass(model, light);
		pointLightPass(model, light);
	}
	glDisable(GL_STENCIL_TEST);

	// handle directional lights
	Transform identity;
	dirlight_pass_shader.use();
	dirlight_pass_shader.setMat4("projection", identity.m);
	dirlight_pass_shader.setMat4("view", identity.m);
	dirlight_pass_shader.setMat4("model", identity.m);

	dirlight_pass_shader.setVec3("viewPos", Vector3f(camera->position));
	dirlight_pass_shader.setVec2("screenSize", width, height);
	dirlight_pass_shader.setFloat("material.shininess", 64.0f);
	for (auto& light : dir_lights) {
		dirLightPass(model, light);
	}
}

void DeferredRenderer::finalPass() {
	gbuffer.bindFinalPass();
	glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_LINEAR);
}

void DeferredRenderer::pointLightPass(GLModel& model, GLPointLight& light) {
	pointlight_pass_shader.use();
	gbuffer.bindLightPass();

	glStencilFunc(GL_NOTEQUAL, 0, 0xFF);

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);

	auto diffuse = light.color * light.intensity;
	pointlight_pass_shader.setBool("light.valid", light.active && light.cast_shadow);
	pointlight_pass_shader.setVec3("light.position", Vector3f(light.position));
	pointlight_pass_shader.setVec3("light.ambient", diffuse * 0.1f);
	pointlight_pass_shader.setVec3("light.diffuse", diffuse);
	pointlight_pass_shader.setVec3("light.specular", 1.0f, 1.0f, 1.0f);
	pointlight_pass_shader.setFloat("light.constant", light.attenuation.constant);
	pointlight_pass_shader.setFloat("light.linear", light.attenuation.linear);
	pointlight_pass_shader.setFloat("light.quadratic", light.attenuation.quadratic);
	pointlight_pass_shader.setFloat("light.far_plane", light.far_plane);

	light.bindShadowMap(pointlight_pass_shader, "light.shadow_map", gbuffer.num_textures);

	// draw bounding sphere of light
	float r = getBoundingSphere(light);
	Transform bs_model = translate(Vector3f(light.position)) * scale(r, r, r);
	pointlight_pass_shader.setMat4("model", bs_model.m);
	bounding_sphere.draw(pointlight_pass_shader, gbuffer.num_textures + 1);

	glCullFace(GL_BACK);
	glDisable(GL_BLEND);
}

void DeferredRenderer::stencilPass(GLModel& model, GLPointLight& light) {
	stencil_pass_shader.use();
	gbuffer.bindStencilPass();

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glClear(GL_STENCIL_BUFFER_BIT);

	glStencilFunc(GL_ALWAYS, 0, 0);	///< always succeeds
	glStencilOpSeparate(GL_BACK, GL_KEEP, GL_INCR_WRAP, GL_KEEP);
	glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_DECR_WRAP, GL_KEEP);

	float r = getBoundingSphere(light);
	Transform bs_model = translate(Vector3f(light.position)) * scale(r, r, r);
	stencil_pass_shader.setMat4("model", bs_model.m);
	bounding_sphere.draw(stencil_pass_shader, 0);
}

void DeferredRenderer::dirLightPass(GLModel& model, GLDirectionalLight& light) {
	dirlight_pass_shader.use();
	gbuffer.bindLightPass();

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);	///< quad is facing the wrong way, so we do this

	auto diffuse = light.color * light.intensity;
	dirlight_pass_shader.setBool("light.valid", light.active && light.cast_shadow);
	dirlight_pass_shader.setVec3("light.direction", light.direction);
	dirlight_pass_shader.setVec3("light.ambient", diffuse * 0.1f);
	dirlight_pass_shader.setVec3("light.diffuse", diffuse);
	dirlight_pass_shader.setVec3("light.specular", 1.0f, 1.0f, 1.0f);
	dirlight_pass_shader.setMat4("light.light_space_mat", light.light_space_mat.m);

	light.bindShadowMap(dirlight_pass_shader, "light.shadow_map", gbuffer.num_textures);
	bounding_quad.draw(dirlight_pass_shader, gbuffer.num_textures + 1);

	glCullFace(GL_BACK);
	glDisable(GL_BLEND);
}

float DeferredRenderer::getBoundingSphere(GLPointLight& light) {
	float max_dim = fmax(fmax(light.color.x, light.color.y), light.color.z);

	float result = (-light.attenuation.linear + sqrtf(light.attenuation.linear * light.attenuation.linear
														- 4 * light.attenuation.quadratic * (light.attenuation.constant
														- 256 * max_dim * light.intensity)))
			/ (2 * light.attenuation.quadratic);

	return result;
}

void DeferredRenderer::generateShadowMaps(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights) {
	// Render depth map
	float near_plane = -5.0f, far_plane = 5.0f;

	depth_shader.use();
	// loop through lights and generate shadow maps
	glViewport(0, 0, shadow_res, shadow_res);
	glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
	for (auto& light : dir_lights) {
		if (light.should_update) {
			depth_shader.setMat4("model", model_mat.m);
			light.generateShadowMap(depth_shader, near_plane, far_plane);
			model.draw(depth_shader, 2); ///< change tex_offset possibly
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}
	}

	// Render depth cube map
	float near = 0.01f, far = 25.0f;

	depth_cube_shader.use();
	// loop through lights and generate shadow maps for point lights
	glViewport(0, 0, shadow_res, shadow_res);
	glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
	for (auto& light : point_lights) {
		if (light.should_update) {
			depth_cube_shader.setMat4("model", model_mat.m);
			light.generateShadowMap(depth_cube_shader, near, far);
			model.draw(depth_cube_shader, 2);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}
	}
}

DeferredRenderer* DeferredRenderer::getSingleton() {
	static DeferredRenderer singleton;
	return &singleton;
}

} // namespace civet