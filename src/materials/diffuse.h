#ifndef CIVET_DIFFUSE_H
#define CIVET_DIFFUSE_H

#include <core/material.h>

namespace civet {

class DiffuseMaterial : public Material {
public:
	DiffuseMaterial(const std::shared_ptr<Texture<Spectrum>>& Kd_, const std::shared_ptr<Texture<float>>& roughness_,
			const std::shared_ptr<Texture<float>>& bump_map_) :
			Kd(Kd_), roughness(roughness_), bump_map(bump_map_) {}

	void computeScatteringFunctions(SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode, bool allow_multiple_lobes) const override;

private:
	std::shared_ptr<Texture<Spectrum>> Kd;
	std::shared_ptr<Texture<float>> roughness, bump_map;
};

} // namespace civet

#endif // CIVET_DIFFUSE_H
