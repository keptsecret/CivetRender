#ifndef CIVET_TEXTURE_SCALE_H
#define CIVET_TEXTURE_SCALE_H

#include <core/civet.h>
#include <core/texture.h>

namespace civet {

template <typename T1, typename T2>
class ScaleTexture : public Texture<T2> {
public:
	ScaleTexture(const std::shared_ptr<Texture<T1>>& t1, const std::shared_ptr<Texture<T2>>& t2) :
			tex1(t1), tex2(t2) {}
	T2 evaluate(const SurfaceInteraction& si, int channel = 0) const {
		return tex1.evaluate(si) * tex2->evaluate(si);
	}

private:
	std::shared_ptr<Texture<T1>> tex1;
	std::shared_ptr<Texture<T2>> tex2;
};

} // namespace civet

#endif // CIVET_TEXTURE_SCALE_H
