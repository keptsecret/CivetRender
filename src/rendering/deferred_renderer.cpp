#include "core/scene.h"
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

	bounding_sphere = GLModel("bounding_sphere");
	bounding_sphere.loadModel("../civet/resources/basic-meshes/sphere.obj");
	bounding_quad = GLModel("bounding_quad");
	bounding_quad.loadModel("../civet/resources/basic-meshes/quad.obj");

	geometry_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_geometry_frag.glsl");
	pointlight_pass_shader = Shader("../civet/src/shaders/deferred_light_pass_vert.glsl", "../civet/src/shaders/deferred_pointlight_pass_frag.glsl");
	dirlight_pass_shader = Shader("../civet/src/shaders/deferred_light_pass_vert.glsl", "../civet/src/shaders/deferred_dirlight_pass_frag.glsl");
	stencil_pass_shader = Shader("../civet/src/shaders/null_vert.glsl", "../civet/src/shaders/null_frag.glsl");
	postprocess_shader = Shader("../civet/src/shaders/deferred_light_pass_vert.glsl", "../civet/src/shaders/postprocess_pass_frag.glsl");

	pointlight_pass_shader.use();
	pointlight_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	pointlight_pass_shader.setInt("AlbedoMap", GBuffer::GBUFFER_TEXTURE_ALBEDO);
	pointlight_pass_shader.setInt("MetallicRoughAOMap", GBuffer::GBUFFER_TEXTURE_METALLICROUGHAO);
	pointlight_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	dirlight_pass_shader.use();
	dirlight_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	dirlight_pass_shader.setInt("AlbedoMap", GBuffer::GBUFFER_TEXTURE_ALBEDO);
	dirlight_pass_shader.setInt("MetallicRoughAOMap", GBuffer::GBUFFER_TEXTURE_METALLICROUGHAO);
	dirlight_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	postprocess_shader.use();
	Transform identity;
	postprocess_shader.use();
	postprocess_shader.setMat4("projection", identity.m);
	postprocess_shader.setMat4("view", identity.m);
	postprocess_shader.setMat4("model", identity.m);

	depth_shader = Shader("../civet/src/shaders/light_depth_vert.glsl", "../civet/src/shaders/light_depth_frag.glsl");
	depth_cube_shader = Shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/light_cube_depth_frag.glsl", "../civet/src/shaders/light_cube_depth_geom.glsl");
}

void DeferredRenderer::draw(Scene& scene) {
	gbuffer.start();

	for (auto model : scene.models) {
		generateShadowMaps(*model, scene.dir_lights, scene.point_lights);

		glViewport(0, 0, width, height);
		glCullFace(GL_BACK);
		geometryPass(*model);
		lightsPass(*model, scene.dir_lights, scene.point_lights);
	}
	postProcessPass(scene);
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
	model.draw(geometry_pass_shader, 0);

	glDepthMask(GL_FALSE);
}

void DeferredRenderer::lightsPass(GLModel& model, std::vector<std::shared_ptr<GLDirectionalLight>>& dir_lights, std::vector<std::shared_ptr<GLPointLight>>& point_lights) {
	glEnable(GL_STENCIL_TEST); ///< render light only if it passes test
	// handle point lights
	stencil_pass_shader.use();
	stencil_pass_shader.setMat4("projection", projection_mat.m);
	stencil_pass_shader.setMat4("view", view_mat.m);

	pointlight_pass_shader.use();
	pointlight_pass_shader.setMat4("projection", projection_mat.m);
	pointlight_pass_shader.setMat4("view", view_mat.m);

	pointlight_pass_shader.setVec3("viewPos", Vector3f(camera->position));
	pointlight_pass_shader.setVec2("screenSize", width, height);
	for (auto light : point_lights) {
		if (light->active) {
			stencilPass(model, *light);
			pointLightPass(model, *light);
		}
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
	for (auto light : dir_lights) {
		if (light->active) {
			dirLightPass(model, *light);
		}
	}
}

void DeferredRenderer::postProcessPass(Scene& scene) {
	scene.skybox->draw(projection_mat, view_mat);

	postprocess_shader.use();
	gbuffer.bindPostProcessPass();

	glDisable(GL_DEPTH_TEST);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT); ///< quad is facing the wrong way, so we do this

	postprocess_shader.setVec2("screenSize", width, height);
	bounding_quad.draw(postprocess_shader, 2);

	glCullFace(GL_BACK);
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

	auto color = light.color * light.power;
	pointlight_pass_shader.setBool("light.valid", light.active && light.cast_shadow);
	pointlight_pass_shader.setVec3("light.position", Vector3f(light.position));
	pointlight_pass_shader.setVec3("light.color", color);
	pointlight_pass_shader.setFloat("light.constant", light.attenuation.constant);
	pointlight_pass_shader.setFloat("light.linear", light.attenuation.linear);
	pointlight_pass_shader.setFloat("light.quadratic", light.attenuation.quadratic);
	pointlight_pass_shader.setFloat("light.radius", light.radius);

	light.bindShadowMap(pointlight_pass_shader, "light.shadow_map", gbuffer.num_textures);

	// draw bounding sphere of light
	float r = getBoundingSphere(light);
	Transform bs_model = translate(Vector3f(light.position)) * scale(r, r, r);
	bounding_sphere.setTransform(bs_model);
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

	glStencilFunc(GL_ALWAYS, 0, 0); ///< always succeeds
	glStencilOpSeparate(GL_BACK, GL_KEEP, GL_INCR_WRAP, GL_KEEP);
	glStencilOpSeparate(GL_FRONT, GL_KEEP, GL_DECR_WRAP, GL_KEEP);

	float r = getBoundingSphere(light);
	Transform bs_model = translate(Vector3f(light.position)) * scale(r, r, r);
	bounding_sphere.setTransform(bs_model);
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
	glCullFace(GL_FRONT); ///< quad is facing the wrong way, so we do this

	auto color = light.color * light.power;
	dirlight_pass_shader.setBool("light.valid", light.active && light.cast_shadow);
	dirlight_pass_shader.setVec3("light.direction", light.direction);
	dirlight_pass_shader.setVec3("light.color", color);
	dirlight_pass_shader.setMat4("light.light_space_mat", light.light_space_mat.m);

	light.bindShadowMap(dirlight_pass_shader, "light.shadow_map", gbuffer.num_textures);
	bounding_quad.draw(dirlight_pass_shader, gbuffer.num_textures + 1);

	glCullFace(GL_BACK);
	glDisable(GL_BLEND);
}

float DeferredRenderer::getBoundingSphere(GLPointLight& light) {
	float max_dim = fmaxf(fmaxf(light.color.x, light.color.y), light.color.z);

	float result = (-light.attenuation.linear + sqrtf(light.attenuation.linear * light.attenuation.linear
														- 4.0f * light.attenuation.quadratic * (light.attenuation.constant
														- 256.0f * max_dim * light.power)))
			/ (2.0f * light.attenuation.quadratic);

	return result;
}

void DeferredRenderer::generateShadowMaps(GLModel& model, std::vector<std::shared_ptr<GLDirectionalLight>>& dir_lights, std::vector<std::shared_ptr<GLPointLight>>& point_lights) {
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glEnable(GL_CULL_FACE);

	// Render depth map
	int max_axis = model.bounds.maximumAxis();
	float expand = std::abs(model.bounds.p_max[max_axis] - model.bounds.p_min[max_axis]) * 0.1; // TODO: reorient frustum to light direction
	Bounds3f frustum = bExpand(model.getWorldBounds(), fmax(expand, 10.f));

	depth_shader.use();
	// loop through lights and generate shadow maps
	glViewport(0, 0, shadow_res, shadow_res);
	glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
	for (auto light : dir_lights) {
		light->generateShadowMap(depth_shader, frustum);
		model.draw(depth_shader, 2); ///< change tex_offset possibly
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	// Render depth cube map
	depth_cube_shader.use();
	// loop through lights and generate shadow maps for point lights
	glViewport(0, 0, shadow_res, shadow_res);
	glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
	for (auto light : point_lights) {
		float near = 0.1f, far = getBoundingSphere(*light);
		light->generateShadowMap(depth_cube_shader, near, far);
		model.draw(depth_cube_shader, 2);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	glDepthMask(GL_FALSE);
}

DeferredRenderer* DeferredRenderer::getSingleton() {
	static DeferredRenderer singleton;
	return &singleton;
}

} // namespace civet