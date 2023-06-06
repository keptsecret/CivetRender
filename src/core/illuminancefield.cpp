#include <core/illuminancefield.h>

#include <core/interaction.h>
#include <core/light.h>
#include <core/primitive.h>
#include <core/scene.h>
#include <utils/parallel.h>
#include <utils/reflection.h>
#include <utils/samplers/random.h>

namespace civet {

IlluminanceField::~IlluminanceField() {
	if (!SG_data_texture_array) {
		glDeleteTextures(1, &SG_data_texture_array);
	}
}

void IlluminanceField::initialize() {
	initializeSGSolver(SGProbe::SG_count, SGDistribution::Spherical);

	const SG* initial_guess = initialGuess();
	SG_sharpness = initial_guess[0].sharpness;
	for (int i = 0; i < SGProbe::SG_count; i++) {
		SG_directions[i] = initial_guess[i].axis;
	}

	probes.resize(probe_grid_size.x * probe_grid_size.y * probe_grid_size.z);
	has_bake_data = false;

	glGenTextures(1, &SG_data_texture_array);
	glBindTexture(GL_TEXTURE_1D_ARRAY, SG_data_texture_array);
	glTexImage2D(GL_TEXTURE_1D_ARRAY, 0, GL_RGB32F, probes.size(), SGProbe::SG_count, 0, GL_RGB, GL_FLOAT, nullptr);
	glTexParameterf(GL_TEXTURE_1D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameterf(GL_TEXTURE_1D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	glBindTexture(GL_TEXTURE_1D_ARRAY, 0);

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
			Spectrum Le = std::exp2(-16.f) * scene.skybox->sampleSky(scene.skybox->atmosphere, ray.d); // expose down value from sky sample
			L += beta * Le;
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
	printf("Baking %zu probes at %llu samples per probe", probes.size(), rays_per_probe);

	RandomSampler sampler(rays_per_probe);
	parallelFor([&](int idx) {
		MemoryArena arena;
		SGProbe* probe = &probes[idx];
		probe->init(rays_per_probe);

		Point3i probe_grid_idx = flatIndexToGridIndex(probe_grid_size, idx);
		Point3f probe_pos{ probe_grid_idx.x * cell_size + corner_position.x,
			probe_grid_idx.y * cell_size + corner_position.y,
			probe_grid_idx.z * cell_size + corner_position.z };

		std::unique_ptr<Sampler> probe_sampler = sampler.clone(idx);
		probe_sampler->startPixel(Point2i());
		uint64_t sample_idx = 0;
		do {
			Ray ray;
			ray.o = probe_pos;
			ray.d = probe->sampleDirection(probe_sampler->get2D());

			Spectrum L(0.f);
			L = pathTrace(ray, scene, *probe_sampler, arena, 5);

			float rgb[3];
			L.toRGB(rgb);

			probe->addSample(ray.d, sample_idx++, Vector3f(rgb[0], rgb[1], rgb[2]));
		} while (probe_sampler->startNextSample());

		probe->bakeResult();
	},
			probes.size(), 1);

	glBindTexture(GL_TEXTURE_1D_ARRAY, SG_data_texture_array);
	for (int i = 0; i < SGProbe::SG_count; i++) {
		for (int j = 0; j < probes.size(); j++) {
			glTexSubImage2D(GL_TEXTURE_1D_ARRAY, 0, j, i, 1, 1, GL_RGB, GL_FLOAT, &probes[j].amplitudes[i]);
		}
	}
	glBindTexture(GL_TEXTURE_1D_ARRAY, 0);

	glCheckError("ERROR::IlluminanceField::bake: OpenGL error code");
	has_bake_data = true;
}

void IlluminanceField::bind(Shader& shader, int tex_offset) {
	shader.setIVec3("probeGridDims", probe_grid_size);
	shader.setVec3("gridCornerCoord", corner_position.x, corner_position.y, corner_position.z);
	shader.setFloat("gridCellSize", cell_size);
	shader.setInt("SGCount", SGProbe::SG_count);
	shader.setFloat("SGSharpness", SG_sharpness);

	for (int i = 0; i < SGProbe::SG_count; i++) {
		shader.setVec3("SGDirections[" + std::to_string(i) + "]", SG_directions[i]);
	}

	glActiveTexture(GL_TEXTURE0 + tex_offset);	///< ensure tex_offset is num textures in gbuffer
	glBindTexture(GL_TEXTURE_1D_ARRAY, SG_data_texture_array);
}

void IlluminanceField::fitGridToBounds(const Bounds3f& bounds) {
	Bounds3f fitbounds = bExpand(bounds, cell_size * 0.5f);
	corner_position = fitbounds.p_min;
	int x_dim = std::ceil(std::abs(fitbounds.p_max.x - fitbounds.p_min.x) / cell_size);
	int y_dim = std::ceil(std::abs(fitbounds.p_max.y - fitbounds.p_min.y) / cell_size);
	int z_dim = std::ceil(std::abs(fitbounds.p_max.z - fitbounds.p_min.z) / cell_size);
	probe_grid_size = Point3i(x_dim, y_dim, z_dim);
}

} // namespace civet