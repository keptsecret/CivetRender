#ifndef CIVET_MATERIAL_H
#define CIVET_MATERIAL_H

#include <core/civet.h>
#include <utils/memory.h>

namespace civet {

enum class TransportMode { Radiance, Importance };

class Material {
public:
	virtual ~Material();

	virtual void computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const = 0;

	static void bump(const std::shared_ptr<Texture<float>>& map, SurfaceInteraction* si);
};

} // namespace civet

#endif // CIVET_MATERIAL_H
