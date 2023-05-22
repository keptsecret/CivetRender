#ifndef CIVET_MEDIUM_H
#define CIVET_MEDIUM_H

#include <core/spectrum.h>

namespace civet {

// TODO: implement ray tracing in medium

class Medium {
public:
	virtual ~Medium() {}
	virtual Spectrum Tr(const Ray& ray, Sampler& sampler) const = 0;
	virtual Spectrum sample(const Ray& ray, Sampler& sampler, MemoryArena& arena, MediumInteraction* mi) const = 0;
};

struct MediumInterface {
	MediumInterface() :
			inside(nullptr), outside(nullptr) {}
	MediumInterface(const Medium* medium) :
			inside(medium), outside(medium) {}
	MediumInterface(const Medium* inside, const Medium* outside) :
			inside(inside), outside(outside) {}

	bool isMediumTransition() const { return inside != outside; }

	const Medium *inside, *outside;
};

} // namespace civet

#endif // CIVET_MEDIUM_H
