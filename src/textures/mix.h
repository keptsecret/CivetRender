#ifndef CIVET_TEXTURE_MIX_H
#define CIVET_TEXTURE_MIX_H

#include <core/civet.h>
#include <core/texture.h>

namespace civet {

template <typename T>
class MixTexture : public Texture<T> {
public:
	MixTexture(const std::shared_ptr<Texture<T>>& t1, const std::shared_ptr<Texture<T>>& t2, const std::shared_ptr<Texture<float>>& amt) :
			tex1(t1), tex2(t2), amount(amt) {}
	T evaluate(const SurfaceInteraction& si, int channel = 0) const {
		T t1 = tex1.evaluate(si), t2 = tex2->evaluate(si);
		float amt = amount->evaluate(si);
		return (1 - amt) * t1 + amt * t2;
	}

private:
	std::shared_ptr<Texture<T>> tex1, tex2;
	std::shared_ptr<Texture<float>> amount;
};

} // namespace civet

#endif // CIVET_TEXTURE_MIX_H
