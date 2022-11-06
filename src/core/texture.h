#ifndef CIVET_TEXTURE_H
#define CIVET_TEXTURE_H

#include <core/civet.h>

namespace civet {

template <typename T>
class Texture {
public:
	virtual T evaluate(const SurfaceInteraction &) const = 0;
	virtual ~Texture() {}
};

} // namespace civet

#endif // CIVET_TEXTURE_H
