#ifndef CIVET_MATERIAL_H
#define CIVET_MATERIAL_H

#include <core/civet.h>

namespace civet {

enum class TransportMode { Radiance, Importance };

class Material {
public:
	virtual ~Material();

	virtual void computeScatteringFunction(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const = 0;
};

} // namespace civet

#endif // CIVET_MATERIAL_H
