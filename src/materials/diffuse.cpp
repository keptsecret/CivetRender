#include <materials/diffuse.h>

#include <utils/reflection.h>
#include <core/interaction.h>
#include <core/texture.h>

namespace civet {

void DiffuseMaterial::computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
	if (bump_map) {
		bump(bump_map, isect);
	}

	isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);
	Spectrum r = Kd->evaluate(*isect).clamp();
	float sig = clamp(roughness->evaluate(*isect), 0, 90);
	if (!r.isBlack()) {
		if (sig == 0) {
			isect->bsdf->add(ARENA_ALLOC(arena, LambertianReflection)(r));
		} else {
			isect->bsdf->add(ARENA_ALLOC(arena, OrenNayar)(r, sig));
		}
	}
}

} // namespace civet