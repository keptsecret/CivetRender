#include <core/skybox.h>

#include <random>
#include <utils/parallel.h>

namespace civet {

const Vector3f Atmosphere::beta_R{3.8e-6f, 13.5e-6f, 33.1e-6f};
const Vector3f Atmosphere::beta_M{21e-6f, 21e-6f, 21e-6f};

bool raySphereIntersect(const Vector3f& o, const Vector3f& d, const float& radius, float& t0, float& t1) {
	// simplified version of the one in spehre.cpp
	float a = d.x * d.x + d.y * d.y + d.z * d.z;
	float b = 2 * (d.x * o.x + d.y * o.y + d.z * o.z);
	float c = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;

	if (!quadratic(a, b, c, t0, t1)) {
		return false;
	}

	if (t0 > t1) {
		std::swap(t0, t1);
	}

	return true;
}

bool Atmosphere::computeIncidentLight(const Vector3f& o, const Vector3f& d, float tmin, float tmax, Vector3f* out, unsigned int num_samples, unsigned int num_light_samples) const {
	float t0, t1;
	if (!raySphereIntersect(o, d, atmosphere_radius, t0, t1) || t1 < 0) {
		return false;
	}

	if (t0 > tmin && t0 > 0) {
		tmin = t0;
	}
	if (t1 < tmax) {
		tmax = t1;
	}

	float segment_length = (tmax - tmin) / num_samples;
	float t = tmin;
	Vector3f sum_R, sum_M;
	float optical_depth_R = 0, optical_depth_M = 0;
	float mu = dot(d, sun_direction); // cosine of the angle between the sun direction and the ray direction
	float phase_R = 3.f / (16.f * Pi) * (1 + mu * mu);
	float g = 0.76f;
	float phase_M = 3.f / (8.f * Pi) * ((1.f - g * g) * (1.f + mu * mu)) / ((2.f + g * g) * pow(1.f + g * g - 2.f * g * mu, 1.5f));

	for (int i = 0; i < num_samples; i++) {
		Vector3f sample_pos = o + (t + segment_length * 0.5f) * d;
		float height = sample_pos.length() - earth_radius;

		// optical depth for light
		float hr = exp(-height / Hr) * segment_length;
		float hm = exp(-height / Hm) * segment_length;
		optical_depth_R += hr;
		optical_depth_M += hm;

		// light optical depth
		float t0light, t1light;
		raySphereIntersect(o, d, atmosphere_radius, t0light, t1light);
		float segment_length_light = t1light / num_light_samples;
		float tlight = 0;
		float optical_depth_light_R = 0, optical_depth_light_M = 0;

		int j;
		for (j = 0; j < num_light_samples; j++) {
			Vector3f sample_pos_light = sample_pos + (tlight + segment_length_light * 0.5f) * sun_direction;
			float height_light = sample_pos_light.length() - earth_radius;
			if (height_light < 0) {
				break;
			}
			optical_depth_light_R += exp(-height_light / Hr) * segment_length_light;
			optical_depth_light_M += exp(-height_light / Hm) * segment_length_light;
			tlight += segment_length_light;
		}
		if (j == num_light_samples) {
			Vector3f tau = beta_R * (optical_depth_R + optical_depth_light_R) +
					beta_M * 1.1f * (optical_depth_M + optical_depth_light_M);
			Vector3f attenuation{ exp(-tau.x), exp(-tau.y), exp(-tau.z) };
			sum_R += attenuation * hr;
			sum_M += attenuation * hm;
		}
		t += segment_length;
	}

	Vector3f result{
		sum_R.x * beta_R.x * phase_R + sum_M.x * beta_M.x * phase_M,
		sum_R.y * beta_R.y * phase_R + sum_M.y * beta_M.y * phase_M,
		sum_R.z * beta_R.z * phase_R + sum_M.z * beta_M.z * phase_M
	};

	*out = result * Skybox::sun_intensity;

	return true;
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
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
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

	Matrix4 view_wot{view.m.m[0][0], view.m.m[0][1], view.m.m[0][2], 0,
		view.m.m[1][0], view.m.m[1][1], view.m.m[1][2], 0,
		view.m.m[2][0], view.m.m[2][1], view.m.m[2][2], 0,
		0, 0, 0, 0};	// strip translate from view matrix

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
	Vector3f sd{ sun_dir };
	sd.y = clamp(sd.y, 0, 1);
	sd = normalize(sd);

	atmosphere.sun_direction = sd;

	float aspect_ratio = 1.f;
	float fov = 90.f;
	int num_samples = parameters.samples_per_pixel;
	for (int i = 0; i < 6; i++) {
		renderSky(i, fov, aspect_ratio, num_samples);
		glTexSubImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, 0, 0, parameters.resolution, parameters.resolution, GL_RGB, GL_FLOAT, skybox_data[i].data());
		glCheckError("ERROR::Skybox::renderSkyboxToTexture: OpenGL error code");
	}
}

void Skybox::renderSky(int face, float fov, float aspect, int spp) {
	Vector3f origin{ 0, atmosphere.earth_radius + 1000, 3000 };
	float angle = std::tan(radians(fov) * 0.5f);
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(0, 1);

	parallelFor2D([&](Point2i xy) {
		unsigned int x = xy.x;
		unsigned int y = xy.y;

		Vector3f* p = skybox_data[face].data();
		p += y * parameters.resolution + x;

		for (unsigned int m = 0; m < spp; m++) {
			for (unsigned int n = 0; n < spp; n++) {
				float rayx = (2 * (x + (m + distribution(generator)) / spp) / float(parameters.resolution) - 1) * aspect * angle;
				float rayy = (1 - (y + (n + distribution(generator)) / spp) / float(parameters.resolution) * 2) * angle;
				Vector3f dir = mapToDirection(rayx, rayy, face);
				float t0, t1, tmax = Infinity;
				if (raySphereIntersect(origin, dir, atmosphere.earth_radius, t0, t1) && t1 > 0) {
					tmax = std::max(0.0f, t0);
				}
				Vector3f color;
				if (atmosphere.computeIncidentLight(origin, dir, 0, tmax, &color, 16, 8)) {
					*p += color;
				} else {
					*p += parameters.ground_color;
				}
			}
		}
		*p *= 1.f / (spp * spp);
	}, Point2i(parameters.resolution, parameters.resolution));
}

Vector3f Skybox::mapToDirection(float x, float y, int s) {
	Vector3f dir;

	// +x, -x, +y, -y, +z, -z
	switch (s)
	{
		case 0:
			dir = normalize(Vector3(1.0f, y, -x));
			break;
		case 1:
			dir = normalize(Vector3(-1.0f, y, x));
			break;
		case 2:
			dir = normalize(Vector3(x, 1.0f, -y));
			break;
		case 3:
			dir = normalize(Vector3(x, -1.0f, y));
			break;
		case 4:
			dir = normalize(Vector3(x, y, 1.0f));
			break;
		case 5:
			dir = normalize(Vector3(-x, y, -1.0f));
			break;
	}

	return dir;
}

const float Skybox::vertices[] = {
	-1.0f,  1.0f, -1.0f,
	-1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,
	1.0f,  1.0f, -1.0f,
	-1.0f,  1.0f, -1.0f,

	-1.0f, -1.0f,  1.0f,
	-1.0f, -1.0f, -1.0f,
	-1.0f,  1.0f, -1.0f,
	-1.0f,  1.0f, -1.0f,
	-1.0f,  1.0f,  1.0f,
	-1.0f, -1.0f,  1.0f,

	1.0f, -1.0f, -1.0f,
	1.0f, -1.0f,  1.0f,
	1.0f,  1.0f,  1.0f,
	1.0f,  1.0f,  1.0f,
	1.0f,  1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,

	-1.0f, -1.0f,  1.0f,
	-1.0f,  1.0f,  1.0f,
	1.0f,  1.0f,  1.0f,
	1.0f,  1.0f,  1.0f,
	1.0f, -1.0f,  1.0f,
	-1.0f, -1.0f,  1.0f,

	-1.0f,  1.0f, -1.0f,
	1.0f,  1.0f, -1.0f,
	1.0f,  1.0f,  1.0f,
	1.0f,  1.0f,  1.0f,
	-1.0f,  1.0f,  1.0f,
	-1.0f,  1.0f, -1.0f,

	-1.0f, -1.0f, -1.0f,
	-1.0f, -1.0f,  1.0f,
	1.0f, -1.0f, -1.0f,
	1.0f, -1.0f, -1.0f,
	-1.0f, -1.0f,  1.0f,
	1.0f, -1.0f,  1.0f
};

} // namespace civet