#include <core/light.h>

#include <core/engine.h>

namespace civet {

bool VisibilityTester::unoccluded(const Scene& scene) const {
	return !scene.intersectP(p0.spawnRayTo(p1));
}

void GLDirectionalLight::generateShadowMap(Shader& shader) {
	Engine* engine = Engine::getSingleton();
	glViewport(0, 0, resolution, resolution);
	if (use_cascaded_shadows) {
		light_space_mat = getLightSpaceMatrices();
		glBindBuffer(GL_UNIFORM_BUFFER, UBO);
		for (int i = 0; i < light_space_mat.size(); i++) {
			Matrix4 m_T = transpose(light_space_mat[i].m);
			glBufferSubData(GL_UNIFORM_BUFFER, i * sizeof(Matrix4), sizeof(Matrix4), m_T.m);
		}
		glBindBuffer(GL_UNIFORM_BUFFER, 0);
	} else {
		light_space_mat[0] = getLightSpaceMatrix(engine->view_camera.near_plane, engine->view_camera.far_plane);
		shader.setMat4("lightSpaceMatrix", light_space_mat[0].m);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glClear(GL_DEPTH_BUFFER_BIT);
	glCheckError("ERROR::GLDirectionalLight::bindShadowMap: OpenGL error code");
	///< maybe change to draw scene in here immediately at some point (pass in scene object)
	///< currently needs to call draw and unbind framebuffer
}

void GLDirectionalLight::bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) {
	auto final_color = color * power;
	shader.setBool("light.valid", active && cast_shadow);
	shader.setVec3("light.direction", direction);
	shader.setVec3("light.color", final_color);
	shader.setBool("light.use_cascaded_shadows", use_cascaded_shadows);
	if (use_cascaded_shadows) {
		shader.setFloat("light.far_plane", Engine::getSingleton()->view_camera.far_plane);
		for (int i = 0; i < cascade_levels.size(); i++) {
			shader.setFloat(("light.cascade_distances[" + std::to_string(i) + "]"), cascade_levels[i]);
		}
		shader.setInt("light.shadow_cascades", tex_offset);
		shader.setInt("light.cascade_count", cascade_levels.size());
		glBindBufferBase(GL_UNIFORM_BUFFER, 1, UBO);
	} else {
		shader.setMat4("light.light_space_mat", light_space_mat[0].m);
		shader.setInt("light.shadow_map", tex_offset);
	}

	glActiveTexture(GL_TEXTURE0 + tex_offset);
	if (use_cascaded_shadows) {
		glBindTexture(GL_TEXTURE_2D_ARRAY, shadow_map);
	} else {
		glBindTexture(GL_TEXTURE_2D, shadow_map);
	}
	glCheckError("ERROR::GLDirectionalLight::bindShadowMap: OpenGL error code");
}

void GLDirectionalLight::init() {
	Engine* engine = Engine::getSingleton();
	float far_plane = engine->view_camera.far_plane;
	cascade_levels = std::vector<float>{ far_plane / 50.0f, far_plane / 25.0f, far_plane / 10.0f };

	glGenFramebuffers(1, &FBO);

	glGenTextures(1, &shadow_map);

	if (use_cascaded_shadows) {
		glGenBuffers(1, &UBO);
		glBindBuffer(GL_UNIFORM_BUFFER, UBO);
		glBufferData(GL_UNIFORM_BUFFER, sizeof(Matrix4) * 4, nullptr, GL_STATIC_DRAW);
		glBindBufferBase(GL_UNIFORM_BUFFER, 1, UBO);
		glBindBuffer(GL_UNIFORM_BUFFER, 0);

		glBindTexture(GL_TEXTURE_2D_ARRAY, shadow_map);
		glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT32F, resolution, resolution,
				cascade_levels.size() + 1, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
		glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		float border_color[] = { 1.0f, 1.0f, 1.0f, 1.0f };
		glTexParameterfv(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BORDER_COLOR, border_color);
		glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, FBO);
		glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, shadow_map, 0);
	} else {
		glBindTexture(GL_TEXTURE_2D, shadow_map);
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
	}

	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);

	int status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE) {
		std::cout << "ERROR::GLDirectionalLight::init: Framebuffer is not complete!";
		throw 0;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glCheckError("ERROR::GLDirectionalLight::init: OpenGL error code");
}

std::vector<Point3f> GLDirectionalLight::getFrustumCornersInWorldSpace(const Transform& projection, const Transform& view) {
	Transform inv_vp = inverse(projection * view);
	std::vector<Point3f> corners;
	for (int x = 0; x < 2; x++) {
		for (int y = 0; y < 2; y++) {
			for (int z = 0; z < 2; z++) {
				Point3f c = inv_vp(Point3f{ 2.0f * x - 1.0f, 2.0f * y - 1.0f, 2.0f * z - 1.0f });
				corners.push_back(c);
			}
		}
	}
	return corners;
}

Transform GLDirectionalLight::getLightSpaceMatrix(const float near_plane, const float far_plane) {
	Engine* engine = Engine::getSingleton();
	const Transform projection = perspective(engine->view_camera.zoom, engine->width / engine->height, near_plane, far_plane);
	const std::vector<Point3f> corners = getFrustumCornersInWorldSpace(projection, engine->view_camera.getViewTransform());

	Point3f center;
	for (const auto& c : corners) {
		center += c;
	}
	center /= corners.size();

	const Transform light_view = lookAtRH(center + direction, center, Vector3f{ 0, 1, 0 });
	Bounds3f frustum{ center, center };
	for (const auto& c : corners) {
		Point3f trf = light_view(c);
		frustum = bUnion(frustum, trf);
	}

	if (frustum.p_min.z < 0) {
		frustum.p_min.z *= frustum_fitting_factor;
	} else {
		frustum.p_min.z /= frustum_fitting_factor;
	}
	if (frustum.p_max.z < 0) {
		frustum.p_max.z /= frustum_fitting_factor;
	} else {
		frustum.p_max.z *= frustum_fitting_factor;
	}

	const Transform light_projection = orthographic(frustum.p_min.x, frustum.p_max.x, frustum.p_min.y, frustum.p_max.y, frustum.p_min.z, frustum.p_max.z);
	return light_projection * light_view;
}

std::vector<Transform> GLDirectionalLight::getLightSpaceMatrices() {
	auto camera = &(Engine::getSingleton()->view_camera);

	std::vector<Transform> mats;
	for (int i = 0; i < cascade_levels.size() + 1; i++) {
		if (i == 0) {
			mats.push_back(getLightSpaceMatrix(camera->near_plane, cascade_levels[i]));
		} else if (i < cascade_levels.size()) {
			mats.push_back(getLightSpaceMatrix(cascade_levels[i - 1], cascade_levels[i]));
		} else {
			mats.push_back(getLightSpaceMatrix(cascade_levels[i - 1], camera->far_plane));
		}
	}
	return mats;
}

void GLPointLight::generateShadowMap(Shader& shader, float near_plane, float far_plane) {
	radius = far_plane;
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
	shader.setVec3("lightPos", position.x, position.y, position.z);
	glViewport(0, 0, resolution, resolution);
	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glClear(GL_DEPTH_BUFFER_BIT);
	///< same for this as above
	glCheckError("ERROR::GLPointLight::generateShadowMap: OpenGL error code");
}

void GLPointLight::bindShadowMap(Shader& shader, const std::string& name, unsigned int tex_offset) {
	auto final_color = color * power;
	shader.setBool("light.valid", active && cast_shadow);
	shader.setVec3("light.position", Vector3f(position));
	shader.setVec3("light.color", final_color);
	shader.setFloat("light.constant", attenuation.constant);
	shader.setFloat("light.linear", attenuation.linear);
	shader.setFloat("light.quadratic", attenuation.quadratic);
	shader.setFloat("light.radius", radius);

	shader.setInt("light.shadow_map", tex_offset);
	glActiveTexture(GL_TEXTURE0 + tex_offset);
	glBindTexture(GL_TEXTURE_CUBE_MAP, shadow_map);
	glCheckError("ERROR::GLPointLight::bindShadowMap: OpenGL error code");
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

	int status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (status != GL_FRAMEBUFFER_COMPLETE) {
		std::cout << "ERROR::GLDirectionalLight::init: Framebuffer is not complete!";
		throw 0;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glCheckError("ERROR::GLPointLight::init: OpenGL error code");
}

} // namespace civet