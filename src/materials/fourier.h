#ifndef CIVET_FOURIER_H
#define CIVET_FOURIER_H

#include <utils/reflection.h>
#include <core/material.h>

namespace civet {

class FourierMaterial : public Material {
public:
	FourierMaterial(const std::string& filename, const std::shared_ptr<Texture<float>>& bump_map_);

	void computeScatteringFunctions(SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode, bool allow_multiple_lobes) const override;

private:
	FourierBSDFTable bsdf_table;
	std::shared_ptr<Texture<float>> bump_map;
};

} // namespace civet

#endif // CIVET_FOURIER_H
