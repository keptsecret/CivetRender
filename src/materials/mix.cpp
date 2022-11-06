#include <materials/mix.h>

#include <core/interaction.h>
#include <core/spectrum.h>
#include <utils/reflection.h>
#include <core/texture.h>

namespace civet {

void MixMaterial::computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
	Spectrum s1 = factor->evaluate(*isect).clamp();
	Spectrum s2 = (Spectrum(1.0f) - s1).clamp();
	SurfaceInteraction si = *isect;
	mat1->computeScatteringFunctions(isect, arena, mode, allow_multiple_lobes);
	mat2->computeScatteringFunctions(&si, arena, mode, allow_multiple_lobes);

	int n1 = isect->bsdf->numComponents();
	int n2 = si.bsdf->numComponents();
	for (int i = 0; i < n1; i++) {
		isect->bsdf->bxdfs[i] = ARENA_ALLOC(arena, ScaledBxDF)(isect->bsdf->bxdfs[i], s1);
	}
	for (int i = 0; i < n2; i++) {
		isect->bsdf->add(ARENA_ALLOC(arena, ScaledBxDF)(si.bsdf->bxdfs[i], s2));
	}
}

} // namespace civet