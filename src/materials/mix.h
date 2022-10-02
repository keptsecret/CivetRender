#ifndef CIVET_MIX_H
#define CIVET_MIX_H

#include <core/material.h>

namespace civet {

class MixMaterial : public Material {
public:
	MixMaterial(const std::shared_ptr<Material>& m1, const std::shared_ptr<Material>& m2, const std::shared_ptr<Texture<Spectrum>>& scale) :
			mat1(m1), mat2(m2), factor(scale) {}

	void computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const override;

private:
	std::shared_ptr<Material> mat1, mat2;
	std::shared_ptr<Texture<Spectrum>> factor;
};

} // namespace civet

#endif // CIVET_MIX_H
