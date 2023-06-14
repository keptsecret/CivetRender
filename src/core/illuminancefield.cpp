#include <core/illuminancefield.h>

#include <core/engine.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <core/interaction.h>
#include <core/light.h>
#include <core/primitive.h>
#include <core/scene.h>
#include <thirdparty/stb/stb_image_write.h>
#include <utils/parallel.h>
#include <utils/reflection.h>
#include <utils/samplers/random.h>

namespace civet {

IlluminanceField::~IlluminanceField() {
	if (!radiance_texture_array) {
		glDeleteTextures(1, &radiance_texture_array);
	}

	if (!distance_texture_array) {
		glDeleteTextures(1, &distance_texture_array);
	}

	if (!sphere_samples_UBO) {
		glDeleteBuffers(1, &sphere_samples_UBO);
	}
}

Vector3f sampleSphereDirection() {
	Vector3f v;
	RNG rng;
	do {
		v.x = rng.uniformFloat() * 2.f - 1.f;
		v.y = rng.uniformFloat() * 2.f - 1.f;
		v.z = rng.uniformFloat() * 2.f - 1.f;
	} while (v.lengthSquared() >= 1.0 || v.lengthSquared() == 0);

	return normalize(v);
}

std::vector<Vector3f> generateDirectionsInSphere(int num_samples) {
	std::vector<Vector3f> directions(num_samples);
	for (int i = 0; i < num_samples; i++) {
		directions[i] = sampleSphereDirection();
	}

	return directions;
}

void IlluminanceField::initialize() {
	glGenFramebuffers(1, &FBO);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);

	radiance_cubemap_data.resize(6);
	for (int i = 0; i < 6; i++) {
		radiance_cubemap_data[i].resize(cubemap_resolution * cubemap_resolution);
	}

	distance_cubemap_data.resize(6);
	for (int i = 0; i < 6; i++) {
		distance_cubemap_data[i].resize(cubemap_resolution * cubemap_resolution);
	}

	probe_radiance.resize(probe_grid_size.x * probe_grid_size.y * probe_grid_size.z);
	for (int i = 0; i < probe_radiance.size(); i++) {
		probe_radiance[i].resize(octahedral_resolution * octahedral_resolution);
	}

	probe_distance.resize(probe_grid_size.x * probe_grid_size.y * probe_grid_size.z);
	for (int i = 0; i < probe_distance.size(); i++) {
		probe_distance[i].resize(octahedral_resolution * octahedral_resolution);
	}
	has_bake_data = false;

	glGenTextures(1, &radiance_cubemap);
	glBindTexture(GL_TEXTURE_CUBE_MAP, radiance_cubemap);
	for (unsigned int i = 0; i < 6; i++) {
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB32F, cubemap_resolution, cubemap_resolution, 0, GL_RGB, GL_FLOAT, nullptr);
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &distance_cubemap);
	glBindTexture(GL_TEXTURE_CUBE_MAP, distance_cubemap);
	for (unsigned int i = 0; i < 6; i++) {
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RG32F, cubemap_resolution, cubemap_resolution, 0, GL_RG, GL_FLOAT, nullptr);
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	glGenTextures(1, &radiance_texture_array);
	glBindTexture(GL_TEXTURE_2D_ARRAY, radiance_texture_array);
	glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGB32F, octahedral_resolution, octahedral_resolution, probe_radiance.size(), 0, GL_RGB, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glGenTextures(1, &distance_texture_array);
	glBindTexture(GL_TEXTURE_2D_ARRAY, distance_texture_array);
	glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RG32F, octahedral_resolution, octahedral_resolution, probe_distance.size(), 0, GL_RG, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glGenTextures(1, &irradiance_texture_array);
	glBindTexture(GL_TEXTURE_2D_ARRAY, irradiance_texture_array);
	glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGB32F, octahedral_resolution, octahedral_resolution, probe_radiance.size(), 0, GL_RGB, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glGenTextures(1, &filtered_distance_texture_array);
	glBindTexture(GL_TEXTURE_2D_ARRAY, filtered_distance_texture_array);
	glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RG32F, octahedral_resolution, octahedral_resolution, probe_distance.size(), 0, GL_RG, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

	int max_samples = std::max(num_irradiance_samples, num_distance_samples);
	std::vector<Vector3f> samples = generateDirectionsInSphere(max_samples);

	glGenBuffers(1, &sphere_samples_UBO);
	glBindBuffer(GL_UNIFORM_BUFFER, sphere_samples_UBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(Vector3f) * max_samples, samples.data(), GL_STATIC_DRAW);
	glBindBufferBase(GL_UNIFORM_BUFFER, 2, sphere_samples_UBO);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);

	irradiance_shader = Shader("../civet/src/shaders/screen_space_vert.glsl", "../civet/src/shaders/irradiance_frag.glsl");

	unsigned int sampledirs_index = glGetUniformBlockIndex(irradiance_shader.ID, "SphereDirectionSamples");
	glUniformBlockBinding(irradiance_shader.ID, sampledirs_index, 2);

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glCheckError("ERROR::IlluminanceField::initialize: OpenGL error code");
}

/*
 * Reduction of the PathIntegrator class from pbrt for tracing light probes
 */
Spectrum estimateDirect(const Interaction& it, const Point2f& u_scattering, const Light& light, const Point2f& u_light,
		const Scene& scene, Sampler& sampler, MemoryArena& arena) {
	Spectrum Ld(0.f);

	Vector3f wi;
	float light_pdf = 0, scattering_pdf = 0;
	VisibilityTester visibility;
	Spectrum Li = light.sample_Li(it, u_light, &wi, &light_pdf, &visibility);
	if (light_pdf > 0 && !Li.isBlack()) {
		Spectrum f;
		if (it.isSurfaceInteraction()) {
			const SurfaceInteraction& isect = (const SurfaceInteraction&)it;
			f = isect.bsdf->f(isect.wo, wi, BSDF_ALL) * absDot(wi, isect.shading.n);
			scattering_pdf = isect.bsdf->pdf(isect.wo, wi, BSDF_ALL);
		} else {
			// TODO: phase function
		}

		if (!f.isBlack()) {
			if (!visibility.unoccluded(scene)) {
				Li = Spectrum(0.f);
			}
			if (!Li.isBlack()) {
				if (isDeltaLight(light.flags)) {
					Ld += f * Li / light_pdf;
				} else {
					float weight = powerHeuristic(1, light_pdf, 1, scattering_pdf);
					Ld += f * Li * weight / light_pdf;
				}
			}
		}
	}

	if (!isDeltaLight(light.flags)) {
		Spectrum f;
		bool sampled_specular = false;
		if (it.isSurfaceInteraction()) {
			BxDFType sampled_type;
			const SurfaceInteraction& isect = (const SurfaceInteraction&)it;
			f = isect.bsdf->sample_f(isect.wo, &wi, u_scattering, &scattering_pdf, BSDF_ALL, &sampled_type);
			f *= absDot(wi, isect.shading.n);
			sampled_specular = sampled_type & BSDF_SPECULAR;
		} else {
			// TODO: medium interaction
		}
		if (!f.isBlack() && scattering_pdf > 0) {
			float weight = 1;
			if (!sampled_specular) {
				light_pdf = light.pdf_Li(it, wi);
				if (light_pdf == 0) {
					return Ld;
				}
				weight = powerHeuristic(1, light_pdf, 1, scattering_pdf);
			}

			SurfaceInteraction light_isect;
			Ray ray = it.spawnRay(wi);
			Spectrum Tr(1.f);
			bool found_isect = scene.intersect(ray, &light_isect);

			// Spectrum Li(0.f);
			// if (found_isect) {
			//	if (light_isect.primitive->getAreaLight() == &light) {
			//		Li = light_isect.Le(-wi);
			//	}
			// }
		}
	}

	return Ld;
}

Spectrum uniformSampleOneLight(const Interaction& it, const Scene& scene, MemoryArena& arena, Sampler& sampler) {
	int n_lights = scene.lights.size();
	if (n_lights == 0) {
		return Spectrum(0.f);
	}
	int light_num = std::min((int)sampler.get1D() * n_lights, n_lights - 1);
	const std::shared_ptr<Light>& light = scene.lights[light_num];

	Point2f u_light = sampler.get2D();
	Point2f u_scattering = sampler.get2D();
	return (float)n_lights * estimateDirect(it, u_scattering, *light, u_light, scene, sampler, arena);
}

Spectrum pathTrace(const Ray& r, const Scene& scene, Sampler& sampler, MemoryArena& arena, int max_depth) {
	Spectrum L(0.f), beta(1.f);
	Ray ray(r);
	bool specular_bounce = false;
	for (int bounces = 0;; bounces++) {
		SurfaceInteraction isect;
		bool found_isect = scene.intersect(ray, &isect);

		if (bounces == 0 || specular_bounce) {
			// TODO: add emitted light
		}

		if (!found_isect) {
			L += beta * scene.skybox->sampleSky(scene.skybox->atmosphere, ray.d);
			break;
		}
		if (bounces >= max_depth) {
			break;
		}

		isect.computeScatteringFunctions(ray, arena, true);
		if (!isect.bsdf) {
			ray = isect.spawnRay(ray.d);
			bounces--;
			continue;
		}

		L += beta * uniformSampleOneLight(isect, scene, arena, sampler);

		Vector3f wo = -ray.d, wi;
		float pdf;
		BxDFType flags;
		Spectrum f = isect.bsdf->sample_f(wo, &wi, sampler.get2D(), &pdf, BSDF_ALL, &flags);
		if (f.isBlack() || pdf == 0.f) {
			break;
		}
		beta *= f * absDot(wi, isect.shading.n) / pdf;
		specular_bounce = (flags & BSDF_SPECULAR) != 0;
		wi = normalize(wi);
		ray = isect.spawnRay(wi);

		///< possibly subsurface here

		if (bounces > 3) {
			float q = std::max(0.05f, 1 - beta.y());
			if (sampler.get1D() < q) {
				break;
			}
			beta /= 1 - q;
		}
	}

	return L;
}

Point3i flatIndexToGridIndex(const Point3i& grid_dims, int idx) {
	return Point3i(idx % grid_dims.x,
			(idx / grid_dims.x) % grid_dims.y,
			idx / (grid_dims.x * grid_dims.y));
}

float signNotZero(float f) {
	return (f >= 0.f) ? 1.f : -1.f;
}

/** Returns a unit vector. Argument o is an octahedral vector packed via octEncode,
	on the [-1, +1] square
	Matches glsl code
 */
Vector3f octDecode(Point2f o) {
	Vector3f v(o.x, o.y, 1.0 - std::abs(o.x) - std::abs(o.y));
	if (v.z < 0.0) {
		v.x = (1.0 - std::abs(v.y)) * signNotZero(v.x);
		v.y = (1.0 - std::abs(v.x)) * signNotZero(v.y);
	}
	return normalize(v);
}

void IlluminanceField::bake(const Scene& scene) {
	printf("Baking %zu probes\n", probe_radiance.size());

	const int num_total_probes = probe_grid_size.x * probe_grid_size.y * probe_grid_size.z;
	GLModel screenspace_quad("bounding_quad");
	screenspace_quad.loadModel("../civet/resources/basic-meshes/quad.obj");

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);
	glViewport(0, 0, octahedral_resolution, octahedral_resolution);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);

	irradiance_shader.use();
	irradiance_shader.setInt("cubemap", 0);

	for (int idx = 0; idx < num_total_probes; idx++) {
		printf("Baking probe %d\r", idx);

		irradiance_shader.setInt("probeIndex", idx);

		RandomSampler radiance_sampler(rays_per_texel_radiance);
		Point3i probe_grid_idx = flatIndexToGridIndex(probe_grid_size, idx);
		Point3f probe_pos{ probe_grid_idx.x * cell_dim.x + corner_position.x,
			probe_grid_idx.y * cell_dim.y + corner_position.y,
			probe_grid_idx.z * cell_dim.z + corner_position.z };

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_CUBE_MAP, radiance_cubemap);

		// Bake radiance cube map
		for (int face = 0; face < 6; face++) {
			std::vector<Vector3f>& face_data = radiance_cubemap_data[face];
			std::fill(face_data.begin(), face_data.end(), Vector3f());

			parallelFor2D([&](Point2i texel) {
				MemoryArena arena;

				Vector3f* p = face_data.data() + texel.y * cubemap_resolution + texel.x;

				std::unique_ptr<Sampler> probe_sampler = radiance_sampler.clone(idx);
				probe_sampler->startPixel(texel);
				do {
					Point2f offset = probe_sampler->get2D();
					Vector3f sample_dir = mapToDirection(texel.x + offset.x, texel.y + offset.y, face);
					Ray ray(probe_pos, normalize(sample_dir));

					Spectrum L(0.f);
					L = pathTrace(ray, scene, *probe_sampler, arena, 5);

					float rgb[3];
					L.toRGB(rgb);

					*p += Vector3f(rgb[0], rgb[1], rgb[2]);
				} while (probe_sampler->startNextSample());

				*p /= rays_per_texel_radiance;
			},
					Point2i(cubemap_resolution, cubemap_resolution));

			glTexSubImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, 0, 0, 0, cubemap_resolution, cubemap_resolution, GL_RGB, GL_FLOAT, radiance_cubemap_data[face].data());
			glCheckError("ERROR::IlluminanceField::bake radiance cubemap: OpenGL error code");
		}

		irradiance_shader.setInt("numSamples", num_irradiance_samples);
		irradiance_shader.setFloat("lobeSize", irradiance_lobe_size);

		glFramebufferTextureLayer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, irradiance_texture_array, 0, idx);
		screenspace_quad.draw(irradiance_shader, 1);
		glCheckError("ERROR::IlluminanceField::bake irradiance: OpenGL error code");

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_CUBE_MAP, distance_cubemap);

		// Bake distance cube map
		for (int face = 0; face < 6; face++) {
			std::vector<Vector2f>& face_data = distance_cubemap_data[face];
			std::fill(face_data.begin(), face_data.end(), Vector2f());

			parallelFor2D([&](Point2i texel) {
				MemoryArena arena;

				Vector2f* p = face_data.data() + texel.y * cubemap_resolution + texel.x;

				Vector3f sample_dir = mapToDirection(texel.x, texel.y, face);
				Ray ray(probe_pos, normalize(sample_dir));

				SurfaceInteraction isect;
				bool found_isect = scene.intersect(ray, &isect);

				float distance_sq = Infinity;
				if (found_isect) {
					distance_sq = (isect.p - probe_pos).lengthSquared();
				}

				*p = Vector2f(std::sqrt(distance_sq), distance_sq);
			},
					Point2i(cubemap_resolution, cubemap_resolution));

			glTexSubImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + face, 0, 0, 0, cubemap_resolution, cubemap_resolution, GL_RG, GL_FLOAT, distance_cubemap_data[face].data());
			glCheckError("ERROR::IlluminanceField::bake distance cubemap: OpenGL error code");
		}

		irradiance_shader.setInt("numSamples", num_distance_samples);
		irradiance_shader.setFloat("lobeSize", distance_lobe_size);

		glFramebufferTextureLayer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, filtered_distance_texture_array, 0, idx);
		screenspace_quad.draw(irradiance_shader, 1);
		glCheckError("ERROR::IlluminanceField::bake filtered distance: OpenGL error code");
	}
	printf("\n");

	glCullFace(GL_BACK);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	has_bake_data = true;
}

void IlluminanceField::bind(Shader& shader, int tex_offset) {
	shader.setIVec3("probeGridDims", probe_grid_size);
	shader.setVec3("gridCornerCoord", corner_position.x, corner_position.y, corner_position.z);
	shader.setVec3("gridCellSize", cell_dim);
	shader.setInt("numIrradianceSamples", num_irradiance_samples);
	shader.setInt("numDistanceSamples", num_distance_samples);
	shader.setFloat("irradianceLobeSize", irradiance_lobe_size);
	shader.setFloat("distanceLobeSize", distance_lobe_size);

	glActiveTexture(GL_TEXTURE0 + tex_offset); ///< ensure tex_offset is num textures in gbuffer
	glBindTexture(GL_TEXTURE_2D_ARRAY, irradiance_texture_array);
	glActiveTexture(GL_TEXTURE0 + tex_offset + 1);
	glBindTexture(GL_TEXTURE_2D_ARRAY, filtered_distance_texture_array);
	glActiveTexture(GL_TEXTURE0);
}

void IlluminanceField::fitGridToBounds(const Bounds3f& bounds) {
	Bounds3f fitbounds(bounds.p_min - cell_dim * 0.5f, bounds.p_max + cell_dim * 0.5f);
	corner_position = fitbounds.p_min;
	int x_dim = std::ceil(std::abs(fitbounds.p_max.x - fitbounds.p_min.x) / cell_dim.x);
	int y_dim = std::ceil(std::abs(fitbounds.p_max.y - fitbounds.p_min.y) / cell_dim.y);
	int z_dim = std::ceil(std::abs(fitbounds.p_max.z - fitbounds.p_min.z) / cell_dim.z);
	probe_grid_size = Point3i(x_dim, y_dim, z_dim);
}

void IlluminanceField::testPathtracer(const Scene& scene) {
	// has problem with left handed coordinates
	printf("started testing pathtrace...");
	Point2i resolution(64, 64);
	GLCamera camera = Engine::getSingleton()->view_camera;

	Bounds2f screen_window(Point2f(-1, -1), Point2f(1, 1));
	Transform camera_to_screen = perspective(45.f, 1e-2f, 1000.0f);
	Transform screen_to_camera = inverse(camera_to_screen);
	Transform screen_to_raster = scale(resolution.x, resolution.y, 1) *
			scale(1 / (screen_window.p_max.x - screen_window.p_min.x), 1 / (screen_window.p_max.y - screen_window.p_min.y), 1) *
			translate(Vector3f(-screen_window.p_min.x, -screen_window.p_max.y, 0));
	Transform raster_to_screen = inverse(screen_to_raster);
	Transform raster_to_camera = screen_to_camera * raster_to_screen;

	Transform world_to_camera = lookAtLH(camera.position, camera.position + camera.front, camera.up);
	Transform camera_to_world = inverse(world_to_camera);

	int num_samples = 1024;
	RandomSampler sampler(num_samples);
	std::vector<unsigned char> pixels(3 * resolution.x * resolution.y, 0);
	parallelFor2D([&](Point2i xy) {
		MemoryArena arena;

		int seed = xy.y * resolution.x + xy.x;
		std::unique_ptr<Sampler> pixel_sampler = sampler.clone(seed);
		pixel_sampler->startPixel(xy);

		unsigned char* p = pixels.data() + 3 * seed;

		Vector3f col;
		do {
			Point2f st = pixel_sampler->get2D();
			Ray ray(Point3f(0, 0, 0), normalize(Vector3f(raster_to_camera(Point3f(xy.x + st.x, xy.y + st.y, 0)))));
			ray = camera_to_world(ray);

			Spectrum L(0.f);
			L = pathTrace(ray, scene, *pixel_sampler, arena, 5);

			float rgb[3];
			L.toRGB(rgb);

			col.x += rgb[0];
			col.y += rgb[1];
			col.z += rgb[2];
		} while (pixel_sampler->startNextSample());

		col /= (float)num_samples;

		p[0] = int(255.99 * std::sqrt(col.x));
		p[1] = int(255.99 * std::sqrt(col.y));
		p[2] = int(255.99 * std::sqrt(col.z));
	},
			resolution);

	stbi_flip_vertically_on_write(true);
	stbi_write_jpg("testpathtrace.jpg", resolution.x, resolution.y, 3, pixels.data(), 80);
	printf("finished pathtrace\n");
}

Vector3f IlluminanceField::mapToDirection(float x, float y, int s) {
	float u = ((x + 0.5f) / (float)cubemap_resolution) * 2.f - 1.f;
	float v = ((y + 0.5f) / (float)cubemap_resolution) * 2.f - 1.f;
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

} // namespace civet