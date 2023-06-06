#include <core/skybox.h>

#include <core/spectrum.h>
#include <utils/hosek/ArHosekSkyModel.h>
#include <utils/parallel.h>
#include <random>

namespace civet {

float angleBetween(const Vector3f& v1, const Vector3f& v2) {
	return std::acos(std::max(dot(v1, v2), 0.00001f));
}

void Atmosphere::init(Vector3f sun_dir, Vector3f ground_albedo, float turbid) {
	sun_dir.y = clamp(sun_dir.y, 0, 1);
	sun_dir = normalize(sun_dir);
	turbid = clamp(turbid, 1.f, 32.f);

	clear();

	float thetaS = angleBetween(sun_dir, Vector3f(0, 1, 0));
	elevation = 0.5 * Pi - thetaS;
	stateR = arhosek_rgb_skymodelstate_alloc_init(turbid, ground_albedo.x, elevation);
	stateG = arhosek_rgb_skymodelstate_alloc_init(turbid, ground_albedo.x, elevation);
	stateB = arhosek_rgb_skymodelstate_alloc_init(turbid, ground_albedo.x, elevation);

	albedo = ground_albedo;
	sun_direction = sun_dir;
	turbidity = turbid;
}

void Atmosphere::clear() {
	if (stateR != nullptr) {
		arhosekskymodelstate_free(stateR);
	}
	if (stateG != nullptr) {
		arhosekskymodelstate_free(stateG);
	}
	if (stateB != nullptr) {
		arhosekskymodelstate_free(stateB);
	}

	turbidity = 0.f;
	albedo = Vector3f(0, 0, 0);
	elevation = 0.f;
	sun_direction = Vector3f(0, 0, 0);
}

void Skybox::init(const SkyboxParameters& params) {
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);

	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), &vertices, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	glCheckError("ERROR::Skybox::init: OpenGL error code");

	update(params);

	shader = Shader("../civet/src/shaders/skybox_vert.glsl", "../civet/src/shaders/skybox_frag.glsl");
}

void Skybox::update(const SkyboxParameters& params) {
	parameters = params;
	editing_params = params;

	glGenTextures(1, &cubemap);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap);
	for (unsigned int i = 0; i < 6; i++) {
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB32F, params.resolution, params.resolution, 0, GL_RGB, GL_FLOAT, nullptr);
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	skybox_data.resize(6);

	for (int i = 0; i < 6; i++) {
		skybox_data[i].resize(params.resolution * params.resolution);
	}

	renderSkyboxToTexture(params.sun_direction);
	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
	glCheckError("ERROR::Skybox::update: OpenGL error code");
}

void Skybox::draw(const Transform& projection, const Transform& view) {
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	Matrix4 view_wot{ view.m.m[0][0], view.m.m[0][1], view.m.m[0][2], 0,
		view.m.m[1][0], view.m.m[1][1], view.m.m[1][2], 0,
		view.m.m[2][0], view.m.m[2][1], view.m.m[2][2], 0,
		0, 0, 0, 0 }; // strip translate from view matrix

	shader.use();
	shader.setInt("skytexture", 0);
	shader.setMat4("view", view_wot);
	shader.setMat4("projection", projection.m);

	glBindVertexArray(VAO);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_CUBE_MAP, cubemap);
	glDrawArrays(GL_TRIANGLES, 0, 36);
	glBindVertexArray(0);

	glDepthFunc(GL_LESS);
	glDisable(GL_DEPTH_TEST);
	glCheckError("ERROR::Skybox::draw: OpenGL error code");
}

void Skybox::renderSkyboxToTexture(const Vector3f& sun_dir) {
	atmosphere.init(sun_dir, parameters.ground_color, parameters.turbidity);

	for (int i = 0; i < 6; i++) {
		parallelFor2D([&](Point2i xy) {
			unsigned int x = xy.x;
			unsigned int y = xy.y;

			Vector3f* p = skybox_data[i].data();
			p += y * parameters.resolution + x;

			Vector3f dir = mapToDirection(x, y, i);
			Spectrum Li = sampleSky(atmosphere, dir);

			float rgb[3];
			Li.toRGB(rgb);
			*p = Vector3f(rgb[0], rgb[1], rgb[2]);
		},
				Point2i(parameters.resolution, parameters.resolution));

		glTexSubImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, 0, 0, parameters.resolution, parameters.resolution, GL_RGB, GL_FLOAT, skybox_data[i].data());
		glCheckError("ERROR::Skybox::renderSkyboxToTexture: OpenGL error code");
	}
}

Spectrum Skybox::sampleSky(const Atmosphere& cache, const Vector3f& sample_dir) {
	float g = angleBetween(sample_dir, cache.sun_direction);
	float t = angleBetween(sample_dir, Vector3f(0, 1, 0));

	Vector3f radiance;
	radiance.x = float(arhosek_tristim_skymodel_radiance(cache.stateR, t, g, 0));
	radiance.y = float(arhosek_tristim_skymodel_radiance(cache.stateR, t, g, 1));
	radiance.z = float(arhosek_tristim_skymodel_radiance(cache.stateR, t, g, 2));
	radiance *= 683.f;

	float rgb[3] = {radiance.x, radiance.y, radiance.z};
	Spectrum Li = Spectrum::fromRGB(rgb);
	return Li;
}

Vector3f Skybox::mapToDirection(float x, float y, int s) {
	float u = ((x + 0.5f) / (float)parameters.resolution) * 2.f - 1.f;
	float v = ((y + 0.5f) / (float)parameters.resolution) * 2.f - 1.f;
	v *= -1.f;

	Vector3f dir;

	// +x, -x, +y, -y, +z, -z
	switch (s) {
		case 0:
			dir = normalize(Vector3(1.0f, v, -u));
			break;
		case 1:
			dir = normalize(Vector3(-1.0f, v, u));
			break;
		case 2:
			dir = normalize(Vector3(u, 1.0f, -v));
			break;
		case 3:
			dir = normalize(Vector3(u, -1.0f, v));
			break;
		case 4:
			dir = normalize(Vector3(u, v, 1.0f));
			break;
		case 5:
			dir = normalize(Vector3(-u, v, -1.0f));
			break;
	}

	return dir;
}

const float Skybox::vertices[] = {
	-1.0f, 1.0f, -1.0f,
	-1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,
	1.0f, 1.0f, -1.0f,
	-1.0f, 1.0f, -1.0f,

	-1.0f, -1.0f, 1.0f,
	-1.0f, -1.0f, -1.0f,
	-1.0f, 1.0f, -1.0f,
	-1.0f, 1.0f, -1.0f,
	-1.0f, 1.0f, 1.0f,
	-1.0f, -1.0f, 1.0f,

	1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,

	-1.0f, -1.0f, 1.0f,
	-1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	1.0f, -1.0f, 1.0f,
	-1.0f, -1.0f, 1.0f,

	-1.0f, 1.0f, -1.0f,
	1.0f, 1.0f, -1.0f,
	1.0f, 1.0f, 1.0f,
	1.0f, 1.0f, 1.0f,
	-1.0f, 1.0f, 1.0f,
	-1.0f, 1.0f, -1.0f,

	-1.0f, -1.0f, -1.0f,
	-1.0f, -1.0f, 1.0f,
	1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,
	-1.0f, -1.0f, 1.0f,
	1.0f, -1.0f, 1.0f
};

} // namespace civet