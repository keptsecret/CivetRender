#ifndef CIVET_TEXTURE_IMAGE_H
#define CIVET_TEXTURE_IMAGE_H

#include <core/civet.h>
#include <core/mipmap.h>
#include <core/texture.h>
#include <map>

namespace civet {

struct TexInfo {
	TexInfo(const std::string& filename, bool do_trilinear, float max_aniso, ImageWrap wrap_mode, float scale, bool gamma) :
			filename(filename),
			do_trilinear(do_trilinear),
			max_aniso(max_aniso),
			wrap_mode(wrap_mode),
			scale(scale),
			gamma(gamma) {}

	std::string filename;
	bool do_trilinear;
	float max_aniso;
	ImageWrap wrap_mode;
	float scale;
	bool gamma;

	bool operator<(const TexInfo& t) const {
		if (filename != t.filename) {
			return filename < t.filename;
		}
		if (do_trilinear != t.do_trilinear) {
			return do_trilinear < t.do_trilinear;
		}
		if (max_aniso != t.max_aniso) {
			return max_aniso < t.max_aniso;
		}
		if (scale != t.scale) {
			return scale < t.scale;
		}
		if (gamma != t.gamma) {
			return !gamma;
		}
		return wrap_mode < t.wrap_mode;
	}
};

template <typename Tmemory, typename Treturn>
class ImageTexture : public Texture<Treturn> {
public:
	ImageTexture(std::unique_ptr<TextureMapping2D> mapping, const std::string& filename,
			bool do_trilinear, float max_aniso, ImageWrap wrap_mode, float scale, bool gamma);

	Treturn evaluate(const SurfaceInteraction& si, int channel = 0) const {
		Vector2f dstdx, dstdy;
		Point2f st = mapping->map(si, &dstdx, &dstdy);
		Tmemory mem = mipmap->lookup(st, dstdx, dstdy);
		Treturn ret;
		convertOut(mem, &ret, channel);
		return ret;
	}

private:
	static MIPMap<Tmemory>* getTexture(const std::string& filename,
			bool do_trilinear, float max_aniso, ImageWrap wrap_mode, float scale, bool gamma);

	static void convertIn(const RGBSpectrum& from, RGBSpectrum* to, float scale, bool gamma) {
		for (int i = 0; i < RGBSpectrum::n_samples; ++i) {
			(*to)[i] = scale * (gamma ? inverseGammaCorrect(from[i]) : from[i]);
		}
	}

	static void convertIn(const RGBSpectrum& from, float* to, float scale, bool gamma) {
		*to = scale * (gamma ? inverseGammaCorrect(from.y()) : from.y());
	}

	static void convertOut(const RGBSpectrum& from, Spectrum* to) {
		float rgb[3];
		from.toRGB(rgb);
		*to = Spectrum::fromRGB(rgb);
	}

	static void convertOut(const RGBSpectrum& from, Treturn* to, int channel) {
		float rgb[3];
		from.toRGB(rgb);
		RGBSpectrum s = Spectrum::fromRGB(rgb);
		*to = s[channel];
	}

	static void convertOut(float from, float* to) {
		*to = from;
	}

	static void clearCache() {
		textures.erase(textures.begin(), textures.end());
	}

	std::unique_ptr<TextureMapping2D> mapping;
	MIPMap<Tmemory>* mipmap;
	static std::map<TexInfo, std::unique_ptr<MIPMap<Tmemory>>> textures;
};

extern template class ImageTexture<float, float>;
extern template class ImageTexture<RGBSpectrum, float>;
extern template class ImageTexture<RGBSpectrum, Spectrum>;

ImageTexture<float, float>* createImageFloatTexture(const Transform& tex2world, const std::string& filename);
ImageTexture<RGBSpectrum, Spectrum>* createImageSpectrumTexture(const Transform& tex2world, const std::string& filename);

} // namespace civet

#endif // CIVET_TEXTURE_IMAGE_H
