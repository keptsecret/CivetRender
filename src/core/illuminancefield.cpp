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

Spectrum pathTraceDirect(const Ray& ray, const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth);

IlluminanceField::~IlluminanceField() {
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
	initializeSGSolver(SGProbe::SG_count, SGDistribution::Spherical);

	const SG* initial_guess = initialGuess();
	SG_sharpness = initial_guess[0].sharpness;
	for (int i = 0; i < SGProbe::SG_count; i++) {
		SG_directions[i] = initial_guess[i].axis;
	}

	const int num_total_probes = probe_grid_size.x * probe_grid_size.y * probe_grid_size.z;
	probes.resize(probe_grid_size.x * probe_grid_size.y * probe_grid_size.z);

	has_bake_data = false;

	glGenTextures(1, &SG_data_texture_array);
	glBindTexture(GL_TEXTURE_1D_ARRAY, SG_data_texture_array);
	glTexImage2D(GL_TEXTURE_1D_ARRAY, 0, GL_RGB32F, SGProbe::SG_count, num_total_probes, 0, GL_RGB, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_1D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_1D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	glBindTexture(GL_TEXTURE_1D_ARRAY, 0);

	glGenFramebuffers(1, &FBO);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);

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

	glGenTextures(1, &filtered_distance_texture_array);
	glBindTexture(GL_TEXTURE_2D_ARRAY, filtered_distance_texture_array);
	glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RG32F, octahedral_resolution, octahedral_resolution, num_total_probes, 0, GL_RG, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

	int max_samples = num_distance_samples;
	std::vector<Vector3f> samples = generateDirectionsInSphere(max_samples);

	glGenBuffers(1, &sphere_samples_UBO);
	glBindBuffer(GL_UNIFORM_BUFFER, sphere_samples_UBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(Vector3f) * max_samples, samples.data(), GL_STATIC_DRAW);
	glBindBufferBase(GL_UNIFORM_BUFFER, 2, sphere_samples_UBO);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);

	distance_cubemap_shader = Shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/precompute_distance_frag.glsl", "../civet/src/shaders/precompute_distance_geom.glsl");
	filtered_distance_shader = Shader("../civet/src/shaders/screen_space_vert.glsl", "../civet/src/shaders/irradiance_frag.glsl");

	unsigned int sampledirs_index = glGetUniformBlockIndex(filtered_distance_shader.ID, "SphereDirectionSamples");
	glUniformBlockBinding(filtered_distance_shader.ID, sampledirs_index, 2);

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

Spectrum uniformSampleAllLights(const Interaction& it, const Scene& scene, MemoryArena& arena, Sampler& sampler) {
	Spectrum L(0.f);
	for (size_t i = 0; i < scene.lights.size(); i++) {
		const std::shared_ptr<Light>& light = scene.lights[i];
		Point2f u_light = sampler.get2D();
		Point2f u_scattering = sampler.get2D();
		L += estimateDirect(it, u_scattering, *light, u_light, scene, sampler, arena);
	}

	return L;
}

Spectrum specularReflect(const Ray& ray, const SurfaceInteraction& isect, const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) {
	Vector3f wo = isect.wo, wi;
	float pdf;
	BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
	Spectrum f = isect.bsdf->sample_f(wo, &wi, sampler.get2D(), &pdf, type);

	const Normal3f& ns = isect.shading.n;
	if (pdf > 0.f && !f.isBlack() && absDot(wi, ns) != 0.f) {
		Ray rd = isect.spawnRay(wi);
		return f * pathTraceDirect(rd, scene, sampler, arena, depth + 1) * absDot(wi, ns) / pdf;
	} else {
		return Spectrum (0.f);
	}
}

Spectrum pathTraceDirect(const Ray& ray, const Scene& scene, Sampler& sampler, MemoryArena& arena, int depth) {
	Spectrum L(0.f);
	SurfaceInteraction isect;
	if (!scene.intersect(ray, &isect)) {
		return scene.skybox->sampleSky(scene.skybox->atmosphere, ray.d);
	}

	isect.computeScatteringFunctions(ray, arena);
	if (!isect.bsdf) {
		return pathTraceDirect(isect.spawnRay(ray.d), scene, sampler, arena, depth);
	}

	if (scene.lights.size() > 0) {
		L += uniformSampleAllLights(isect, scene, arena, sampler);
	}
	if (depth + 1 < 5) {
		L += specularReflect(ray, isect, scene, sampler, arena, depth);
	}

	return L;
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

void IlluminanceField::bake(const Scene& scene) {
	const int num_total_probes = probe_grid_size.x * probe_grid_size.y * probe_grid_size.z;
	printf("Baking %d probe SGs\n", num_total_probes);

	RandomSampler sampler(rays_per_probe);
	parallelFor([&](int idx) {
		MemoryArena arena;
		SGProbe* probe = &probes[idx];
		probe->init(rays_per_probe);

		Point3i probe_grid_idx = flatIndexToGridIndex(probe_grid_size, idx);
		Point3f probe_pos{ probe_grid_idx.x * cell_dim.x + corner_position.x,
			probe_grid_idx.y * cell_dim.y + corner_position.y,
			probe_grid_idx.z * cell_dim.z + corner_position.z };

		std::unique_ptr<Sampler> probe_sampler = sampler.clone(idx);
		probe_sampler->startPixel(Point2i());
		uint64_t sample_idx = 0;
		do {
			Ray ray(probe_pos, normalize(probe->sampleDirection(probe_sampler->get2D())));

			Spectrum L(0.f);
			L = pathTrace(ray, scene, *probe_sampler, arena, 5);

			float rgb[3];
			L.toRGB(rgb);

			probe->addSample(ray.d, sample_idx++, Vector3f(rgb[0], rgb[1], rgb[2]));
		} while (probe_sampler->startNextSample());

		probe->bakeResult();
	}, num_total_probes, 1);

	glBindTexture(GL_TEXTURE_1D_ARRAY, SG_data_texture_array);
	for (int i = 0; i < probes.size(); i++) {
		glTexSubImage2D(GL_TEXTURE_1D_ARRAY, 0, 0, i, SGProbe::SG_count, 1, GL_RGB, GL_FLOAT, probes[i].amplitudes.data());
	}
	glBindTexture(GL_TEXTURE_1D_ARRAY, 0);

	glCheckError("ERROR::IlluminanceField::bake SGs: OpenGL error code");

	GLModel screenspace_quad("screenspace_quad");
	screenspace_quad.loadModel("../civet/resources/basic-meshes/quad.obj");

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FBO);
	glEnable(GL_CULL_FACE);

	for (int idx = 0; idx < num_total_probes; idx++) {
		printf("Baking probe %d filtered distance\r", idx);

		Point3i probe_grid_idx = flatIndexToGridIndex(probe_grid_size, idx);
		Point3f probe_pos{ probe_grid_idx.x * cell_dim.x + corner_position.x,
			probe_grid_idx.y * cell_dim.y + corner_position.y,
			probe_grid_idx.z * cell_dim.z + corner_position.z };

		// Bake distance cube map
		glViewport(0, 0, cubemap_resolution, cubemap_resolution);
		glEnable(GL_DEPTH_TEST);
		glCullFace(GL_BACK);

		distance_cubemap_shader.use();
		std::vector<Transform> view_transforms;
		view_transforms.push_back(lookAtRH(probe_pos, probe_pos + Vector3f(1, 0, 0), Vector3f(0, -1, 0)));
		view_transforms.push_back(lookAtRH(probe_pos, probe_pos + Vector3f(-1, 0, 0), Vector3f(0, -1, 0)));
		view_transforms.push_back(lookAtRH(probe_pos, probe_pos + Vector3f(0, 1, 0), Vector3f(0, 0, 1)));
		view_transforms.push_back(lookAtRH(probe_pos, probe_pos + Vector3f(0, -1, 0), Vector3f(0, 0, -1)));
		view_transforms.push_back(lookAtRH(probe_pos, probe_pos + Vector3f(0, 0, 1), Vector3f(0, -1, 0)));
		view_transforms.push_back(lookAtRH(probe_pos, probe_pos + Vector3f(0, 0, -1), Vector3f(0, -1, 0)));

		for (unsigned int i = 0; i < 6; i++) {
			distance_cubemap_shader.setMat4("viewMatrices[" + std::to_string(i) + "]", view_transforms[i].m);
		}

		glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, distance_cubemap, 0);

		for (const auto& model : scene.models) {
			model->draw(distance_cubemap_shader, 2);
		}

		// Bake filtered distance
		glViewport(0, 0, octahedral_resolution, octahedral_resolution);
		glDisable(GL_DEPTH_TEST);
		glCullFace(GL_FRONT);

		filtered_distance_shader.use();
		filtered_distance_shader.setInt("cubemap", 0);
		filtered_distance_shader.setInt("probeIndex", idx);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_CUBE_MAP, distance_cubemap);

		filtered_distance_shader.setInt("numSamples", num_distance_samples);
		filtered_distance_shader.setFloat("lobeSize", distance_lobe_size);

		glFramebufferTextureLayer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, filtered_distance_texture_array, 0, idx);
		screenspace_quad.draw(filtered_distance_shader, 1);
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
	shader.setInt("SGCount", SGProbe::SG_count);
	shader.setFloat("SGSharpness", SG_sharpness);
	shader.setFloat("distanceLobeSize", distance_lobe_size);

	for (int i = 0; i < SGProbe::SG_count; i++) {
		shader.setVec3("SGDirections[" + std::to_string(i) + "]", SG_directions[i]);
	}

	glActiveTexture(GL_TEXTURE0 + tex_offset); ///< ensure tex_offset is num textures in gbuffer
	glBindTexture(GL_TEXTURE_1D_ARRAY, SG_data_texture_array);
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

		p[0] = int(255.99 * clamp(std::sqrt(col.x), 0, 1));
		p[1] = int(255.99 * clamp(std::sqrt(col.y), 0, 1));
		p[2] = int(255.99 * clamp(std::sqrt(col.z), 0, 1));
	},
			resolution);

	stbi_flip_vertically_on_write(true);
	stbi_write_jpg("testpathtrace.jpg", resolution.x, resolution.y, 3, pixels.data(), 80);
	printf("finished pathtrace\n");
}

} // namespace civet