#ifndef CIVET_DISNEY_H
#define CIVET_DISNEY_H

#include <core/material.h>

namespace civet {

class DisneyMaterial : public Material {
public:
	DisneyMaterial(const std::shared_ptr<Texture<Spectrum>>& color,
			const std::shared_ptr<Texture<float>>& metallic,
			const std::shared_ptr<Texture<float>>& eta,
			const std::shared_ptr<Texture<float>>& roughness,
			const std::shared_ptr<Texture<float>>& spec_tint,
			const std::shared_ptr<Texture<float>>& anis,
			const std::shared_ptr<Texture<float>>& sheen,
			const std::shared_ptr<Texture<float>>& sheen_tint,
			const std::shared_ptr<Texture<float>>& cc,
			const std::shared_ptr<Texture<float>>& cc_gloss,
			const std::shared_ptr<Texture<float>>& spec_trans,
			const std::shared_ptr<Texture<Spectrum>>& scatter_dist,
			bool t,
			const std::shared_ptr<Texture<float>>& flatness,
			const std::shared_ptr<Texture<float>>& diff_trans,
			const std::shared_ptr<Texture<float>>& bump) :
			color(color),
			metallic(metallic),
			eta(eta),
			roughness(roughness),
			specular_tint(spec_tint),
			anisotropic(anis),
			sheen(sheen),
			sheen_tint(sheen_tint),
			clearcoat(cc),
			clearcoat_gloss(cc_gloss),
			specular_transmission(spec_trans),
			scatter_distance(scatter_dist),
			thin(t),
			flatness(flatness),
			diffuse_transmission(diff_trans),
			bump_map(bump) {}

	void computeScatteringFunctions(SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode, bool allow_multiple_lobes) const override;

private:
	std::shared_ptr<Texture<Spectrum>> color;
	std::shared_ptr<Texture<float>> metallic, eta;
	std::shared_ptr<Texture<float>> roughness, specular_tint, anisotropic, sheen;
	std::shared_ptr<Texture<float>> sheen_tint, clearcoat, clearcoat_gloss;
	std::shared_ptr<Texture<float>> specular_transmission;
	std::shared_ptr<Texture<Spectrum>> scatter_distance;
	bool thin;
	std::shared_ptr<Texture<float>> flatness, diffuse_transmission, bump_map;
};

} // namespace civet

#endif // CIVET_DISNEY_H
