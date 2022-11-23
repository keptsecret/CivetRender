#include <rendering/deferred_renderer.h>

namespace civet {

void DeferredRenderer::init(unsigned int w, unsigned int h) {
	width = w, height = h;
	if (!gbuffer.init(width, height)) {
		printf("ERROR::DeferredRenderer: GBuffer init error");
	}

	bounding_sphere = GLModel("../civet/resources/basic-meshes/sphere.obj");
	bounding_quad = GLModel("../civet/resources/basic-meshes/quad.obj");

	geometry_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_geometry_frag.glsl");
	pointlight_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_pointlight_pass_frag.glsl");
	dirlight_pass_shader = Shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_dirlight_pass_frag.glsl");

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
}

void DeferredRenderer::draw(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights) {
	geometryPass(model);
	lightsPass(model, dir_lights, point_lights);
}

void DeferredRenderer::geometryPass(GLModel& model) {
	geometry_pass_shader.use();
	gbuffer.bindWrite();

	glDepthMask(GL_TRUE); // only update depth buffer on geometry pass
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	geometry_pass_shader.use();
	geometry_pass_shader.setMat4("projection", projection_mat.m);
	geometry_pass_shader.setMat4("view", view_mat.m);
	geometry_pass_shader.setMat4("model", model_mat.m);
	model.draw(geometry_pass_shader, 0);

	glDepthMask(GL_FALSE);
	glDisable(GL_DEPTH_TEST);
}

void DeferredRenderer::lightsPass(GLModel& model, std::vector<GLDirectionalLight>& dir_lights, std::vector<GLPointLight>& point_lights) {
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);

	gbuffer.bindRead();
	glClear(GL_COLOR_BUFFER_BIT);

	// handle point lights
	pointlight_pass_shader.use();
	pointlight_pass_shader.setMat4("projection", projection_mat.m);
	pointlight_pass_shader.setMat4("view", view_mat.m);

	pointlight_pass_shader.setVec3("viewPos", Vector3f(camera->position));
	pointlight_pass_shader.setVec2("screenSize", width, height);
	pointlight_pass_shader.setFloat("material.shininess", 64.0f);
	pointLightsPass(model, point_lights);

	// handle directional lights
	Transform identity;
	dirlight_pass_shader.use();
	dirlight_pass_shader.setMat4("projection", identity.m);
	dirlight_pass_shader.setMat4("view", identity.m);
	dirlight_pass_shader.setMat4("model", identity.m);

	dirlight_pass_shader.setVec3("viewPos", Vector3f(camera->position));
	dirlight_pass_shader.setVec2("screenSize", width, height);
	dirlight_pass_shader.setFloat("material.shininess", 64.0f);
	dirLightsPass(model, dir_lights);
}

void DeferredRenderer::pointLightsPass(GLModel& model, std::vector<GLPointLight>& lights) {
	for (int i = 0; i < lights.size(); i++) {
		auto diffuse = lights[i].color * lights[i].intensity;

		pointlight_pass_shader.setBool("light.valid", lights[i].active && lights[i].cast_shadow);
		pointlight_pass_shader.setVec3("light.ambient", diffuse * 0.1f);
		pointlight_pass_shader.setVec3("light.position", Vector3f(lights[i].position));
		pointlight_pass_shader.setVec3("light.diffuse", diffuse);
		pointlight_pass_shader.setVec3("light.specular", 1.0f, 1.0f, 1.0f);
		pointlight_pass_shader.setFloat("light.constant", lights[i].attenuation.constant);
		pointlight_pass_shader.setFloat("light.linear", lights[i].attenuation.quadratic);
		pointlight_pass_shader.setFloat("light.quadratic", lights[i].attenuation.quadratic);

		pointlight_pass_shader.setFloat("light.far_plane", lights[i].far_plane);

		// draw bounding sphere of light
		float r = getBoundingSphere(lights[i]);
		Transform bs_model = translate(Vector3f(lights[i].position)) * scale(r, r, r); ///< verify correct order of transformation
		pointlight_pass_shader.setMat4("model", bs_model.m);
		bounding_sphere.draw(pointlight_pass_shader, 0);	///< there shouldn't be any shadows?
	}
}

void DeferredRenderer::dirLightsPass(GLModel& model, std::vector<GLDirectionalLight>& lights) {
	for (int i = 0; i < lights.size(); i++) {
		auto diffuse = lights[i].color * lights[i].intensity;

		dirlight_pass_shader.setBool("light.valid", lights[i].active && lights[i].cast_shadow);
		dirlight_pass_shader.setVec3("light.direction", lights[i].direction);
		dirlight_pass_shader.setVec3("light.ambient", diffuse * 0.1f);
		dirlight_pass_shader.setVec3("light.diffuse", diffuse);
		dirlight_pass_shader.setVec3("light.specular", 1.0f, 1.0f, 1.0f);

		dirlight_pass_shader.setMat4("light.light_space_mat", lights[i].light_space_mat.m);

		// TODO: draw screen space quad
		bounding_quad.draw(dirlight_pass_shader, 0);
	}
}

float DeferredRenderer::getBoundingSphere(GLPointLight& light) {
	float max_dim = fmax(fmax(light.color.x, light.color.y), light.color.z);

	float result = (-light.attenuation.linear + sqrtf(light.attenuation.linear * light.attenuation.linear
														- 4 * light.attenuation.quadratic * (light.attenuation.constant
														- 256 * max_dim * light.intensity)))
			/ (2 * light.attenuation.quadratic);

	return result;
}

DeferredRenderer* DeferredRenderer::getSingleton() {
	static DeferredRenderer singleton;
	return &singleton;
}

} // namespace civet