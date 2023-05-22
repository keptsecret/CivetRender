#ifndef CIVET_TEXTURE_CONSTANT_H
#define CIVET_TEXTURE_CONSTANT_H

#include <core/civet.h>
#include <core/texture.h>

namespace civet {

template <typename T>
class ConstantTexture : public Texture<T> {
public:
	ConstantTexture(const T& value) : value(value) {}
	T evaluate(const SurfaceInteraction&) const {
		return value;
	}

private:
	T value;
};

} // namespace civet

#endif // CIVET_TEXTURE_CONSTANT_H
