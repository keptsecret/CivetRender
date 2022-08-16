#include <core/light.h>

namespace civet {

void GLDirectionalLight::generateShadowMap(Shader& shader, float near_plane, float far_plane) {
	Transform projection = orthographic(-10.0f, 10.0f, -10.0f, 10.0f, near_plane, far_plane);
	Transform view = lookAtRH(Point3f(2, -1, -1), Point3f(0, 0, 0), Vector3f(0, 1, 0));
	Transform light_space_tf = projection * view;

	shader.setMat4("lightSpaceMatrix", light_space_tf.m);
	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glClear(GL_DEPTH_BUFFER_BIT);
	///< maybe change to draw scene in here immediately at some point (pass in scene object)
	///< currently needs to call draw and unbind framebuffer
}

void GLDirectionalLight::bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) {
	shader.setInt(name, tex_offset);

	glActiveTexture(GL_TEXTURE0 + tex_offset);
	glBindTexture(GL_TEXTURE_2D, shadow_map);
}

void GLDirectionalLight::init() {
	glGenFramebuffers(1, &FBO);
	glGenTextures(1, &shadow_map);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, resolution, resolution, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	float border_color[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border_color);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadow_map, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void GLPointLight::generateShadowMap(Shader& shader, float near_plane, float far_plane) {
	Transform point_projection = perspective(90.0f, 1.0f, near_plane, far_plane);
	std::vector<Transform> point_transforms;
	point_transforms.push_back(point_projection * lookAtRH(position, position + Vector3f(1, 0, 0), Vector3f(0, -1, 0)));
	point_transforms.push_back(point_projection * lookAtRH(position, position + Vector3f(-1, 0, 0), Vector3f(0, -1, 0)));
	point_transforms.push_back(point_projection * lookAtRH(position, position + Vector3f(0, 1, 0), Vector3f(0, 0, 1)));
	point_transforms.push_back(point_projection * lookAtRH(position, position + Vector3f(0, -1, 0), Vector3f(0, 0, -1)));
	point_transforms.push_back(point_projection * lookAtRH(position, position + Vector3f(0, 0, 1), Vector3f(0, -1, 0)));
	point_transforms.push_back(point_projection * lookAtRH(position, position + Vector3f(0, 0, -1), Vector3f(0, -1, 0)));

	for (unsigned int i = 0; i < 6; i++) {
		shader.setMat4("shadowMatrices[" + std::to_string(i) + "]", point_transforms[i].m);
	}
	shader.setFloat("far_plane", far_plane);
	shader.setVec3("lightPos", Vector3f(position));
	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glClear(GL_DEPTH_BUFFER_BIT);
	///< same for this as above
}

void GLPointLight::bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) {
	shader.setInt(name, tex_offset);

	glActiveTexture(GL_TEXTURE0 + tex_offset);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow_map);
}

void GLPointLight::init() {
	glGenFramebuffers(1, &FBO);
	glGenTextures(1, &shadow_map);

	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow_map);
	for (unsigned int i = 0; i < 6; i++) {
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_DEPTH_COMPONENT, resolution, resolution, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, shadow_map, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

} // namespace civet