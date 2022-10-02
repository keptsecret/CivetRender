#ifndef CIVET_DIELECTRIC_H
#define CIVET_DIELECTRIC_H

#include <core/material.h>

namespace civet {

class DielectricMaterial : public Material {
public:
	DielectricMaterial(const std::shared_ptr<Texture<Spectrum>>& Kd_, const std::shared_ptr<Texture<Spectrum>>& Ks_,
			const std::shared_ptr<Texture<float>>& roughness_, const std::shared_ptr<Texture<float>>& bump_map_, bool remap_roughness_) :
			Kd(Kd_), Ks(Ks_), roughness(roughness_), bump_map(bump_map_), remap_roughness(remap_roughness_) {}

	void computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const override;

private:
	std::shared_ptr<Texture<Spectrum>> Kd, Ks;
	std::shared_ptr<Texture<float>> roughness, bump_map;
	const bool remap_roughness;
};

} // namespace civet

#endif // CIVET_DIELECTRIC_H
