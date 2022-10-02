#include <materials/dielectric.h>

#include <core/texture.h>
#include <utils/microfacet.h>
#include <utils/reflection.h>

namespace civet {

void DielectricMaterial::computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
	if (bump_map) {
		bump(bump_map, isect);
	}

	isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);
	Spectrum kd = Kd->evaluate(*isect).clamp();
	if (!kd.isBlack()) {
		isect->bsdf->add(ARENA_ALLOC(arena, LambertianReflection)(kd));
	}

	Spectrum ks = Ks->evaluate(*isect).clamp();
	if (!ks.isBlack()) {
		Fresnel* fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.0f, 1.5f);
		float rough = roughness->evaluate(*isect);
		if (remap_roughness) {
			rough = TrowbridgeReitzDistribution::roughnessToAlpha(rough);
		}
		MicrofacetDistribution* distrib = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(rough, rough);
		BxDF* specular = ARENA_ALLOC(arena, MicrofacetReflection)(ks, distrib, fresnel);
		isect->bsdf->add(specular);
	}
}

} // namespace civet