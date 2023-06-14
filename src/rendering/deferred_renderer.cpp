#include <rendering/deferred_renderer.h>

#include <core/engine.h>

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
	indirect_pass_shader = Shader("../civet/src/shaders/deferred_light_pass_vert.glsl", "../civet/src/shaders/deferred_indirect_pass_frag.glsl");
	postprocess_shader = Shader("../civet/src/shaders/deferred_light_pass_vert.glsl", "../civet/src/shaders/postprocess_pass_frag.glsl");

	pointlight_pass_shader.use();
	pointlight_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	pointlight_pass_shader.setInt("AlbedoMap", GBuffer::GBUFFER_TEXTURE_ALBEDO);
	pointlight_pass_shader.setInt("AORoughMetallicMap", GBuffer::GBUFFER_TEXTURE_AOROUGHMETALLIC);
	pointlight_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	dirlight_pass_shader.use();
	dirlight_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	dirlight_pass_shader.setInt("AlbedoMap", GBuffer::GBUFFER_TEXTURE_ALBEDO);
	dirlight_pass_shader.setInt("AORoughMetallicMap", GBuffer::GBUFFER_TEXTURE_AOROUGHMETALLIC);
	dirlight_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	indirect_pass_shader.use();
	indirect_pass_shader.setInt("PositionMap", GBuffer::GBUFFER_TEXTURE_POSITION);
	indirect_pass_shader.setInt("AlbedoMap", GBuffer::GBUFFER_TEXTURE_ALBEDO);
	indirect_pass_shader.setInt("AORoughMetallicMap", GBuffer::GBUFFER_TEXTURE_AOROUGHMETALLIC);
	indirect_pass_shader.setInt("NormalMap", GBuffer::GBUFFER_TEXTURE_NORMAL);

	indirect_pass_shader.setInt("radianceOctMap", gbuffer.num_textures);
	indirect_pass_shader.setInt("distanceOctMap", gbuffer.num_textures + 1);

	postprocess_shader.use();
	Transform identity;
	postprocess_shader.use();
	postprocess_shader.setMat4("projection", identity.m);
	postprocess_shader.setMat4("view", identity.m);
	postprocess_shader.setMat4("model", identity.m);

	depth_shader = Shader("../civet/src/shaders/light_depth_vert.glsl", "../civet/src/shaders/light_depth_frag.glsl");
	depth_cascade_shader = Shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/light_depth_frag.glsl", "../civet/src/shaders/light_cascaded_depth_geom.glsl");
	depth_cube_shader = Shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/light_cube_depth_frag.glsl", "../civet/src/shaders/light_cube_depth_geom.glsl");

	// bind light matrices UBO
	unsigned int lightmats_index = glGetUniformBlockIndex(depth_cascade_shader.ID, "LightSpaceMatrices");
	glUniformBlockBinding(depth_cascade_shader.ID, lightmats_index, 1);
	lightmats_index = glGetUniformBlockIndex(dirlight_pass_shader.ID, "LightSpaceMatrices");
	glUniformBlockBinding(dirlight_pass_shader.ID, lightmats_index, 1);
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
	indirectLightingPass(scene);
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

void DeferredRenderer::indirectLightingPass(Scene& scene) {
	if (!scene.probe_grid->hasBakeData()) {
		return;
	}

	indirect_pass_shader.use();
	gbuffer.bindLightPass();

	Transform identity;
	indirect_pass_shader.setMat4("projection", identity.m);
	indirect_pass_shader.setMat4("view", identity.m);
	indirect_pass_shader.setMat4("model", identity.m);
	indirect_pass_shader.setVec3("viewPos", Vector3f(camera->position));
	indirect_pass_shader.setVec2("screenSize", width, height);

	scene.probe_grid->bind(indirect_pass_shader, gbuffer.num_textures);

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT); ///< quad is facing the wrong way, so we do this

	bounding_quad.draw(indirect_pass_shader, gbuffer.num_textures + 2);

	glCullFace(GL_BACK);
	glDisable(GL_BLEND);
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

	light.bindShadowMap(pointlight_pass_shader, "light.shadow_map", gbuffer.num_textures);

	// draw bounding sphere of light
	float r = getBoundingSphere(light);
	bounding_sphere.transform_data.scale_vec = Vector3f{r, r, r};
	bounding_sphere.transform_data.translation = Vector3f{light.position};
	bounding_sphere.transform_data.updateTransform();
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
	bounding_sphere.transform_data.scale_vec = Vector3f{r, r, r};
	bounding_sphere.transform_data.translation = Vector3f{light.position};
	bounding_sphere.transform_data.updateTransform();
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

	dirlight_pass_shader.setMat4("cam_view", view_mat.m);
	light.bindShadowMap(dirlight_pass_shader,  "light.shadow_map", gbuffer.num_textures);
	bounding_quad.draw(dirlight_pass_shader, gbuffer.num_textures + 1);

	glCullFace(GL_BACK);
	glDisable(GL_BLEND);
}

float DeferredRenderer::getBoundingSphere(GLPointLight& light) {
	float max_dim = fmaxf(fmaxf(light.color.x, light.color.y), light.color.z);

	float result = (-light.attenuation.linear + sqrtf(light.attenuation.linear * light.attenuation.linear
														- 4.0f * light.attenuation.quadratic * (light.attenuation.constant
														- 64.0f * max_dim)))	// 256 / 4
			/ (2.0f * light.attenuation.quadratic);

	return result;
}

void DeferredRenderer::generateShadowMaps(GLModel& model, std::vector<std::shared_ptr<GLDirectionalLight>>& dir_lights, std::vector<std::shared_ptr<GLPointLight>>& point_lights) {
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glEnable(GL_CULL_FACE);

	// Render depth map
	// loop through lights and generate shadow maps
	glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
	glEnable(GL_DEPTH_CLAMP);
	for (auto light : dir_lights) {
		if (light->active) {
			if (light->use_cascaded_shadows) {
				depth_cascade_shader.use();
				light->generateShadowMap(depth_cascade_shader);
				model.draw(depth_cascade_shader, 2); ///< change tex_offset possibly
			} else {
				depth_shader.use();
				light->generateShadowMap(depth_shader);
				model.draw(depth_shader, 2); ///< change tex_offset possibly
			}
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}
	}
	glDisable(GL_DEPTH_CLAMP);

	// Render depth cube map
	depth_cube_shader.use();
	// loop through lights and generate shadow maps for point lights
	glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
	for (auto light : point_lights) {
		if (light->active) {
			float near = 0.1f, far = getBoundingSphere(*light);
			light->generateShadowMap(depth_cube_shader, near, far);
			model.draw(depth_cube_shader, 2);
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}
	}

	glDepthMask(GL_FALSE);
}

DeferredRenderer* DeferredRenderer::getSingleton() {
	static DeferredRenderer singleton;
	return &singleton;
}

} // namespace civet