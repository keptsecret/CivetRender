#include <materials/fourier.h>

namespace civet {

FourierMaterial::FourierMaterial(const std::string& filename, const std::shared_ptr<Texture<float>>& bump_map_) :
	bump_map(bump_map_) {
	FourierBSDFTable::read(filename, &bsdf_table);
}

void FourierMaterial::computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
	if (bump_map) {
		bump(bump_map, isect);
	}

	isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);
	isect->bsdf->add(ARENA_ALLOC(arena, FourierBSDF)(bsdf_table, mode));
}

} // namespace civet