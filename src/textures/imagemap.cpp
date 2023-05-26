#include <textures/imagemap.h>

#include <utils/imageio.h>

namespace civet {

template <typename Tmemory, typename Treturn>
ImageTexture<Tmemory, Treturn>::ImageTexture(std::unique_ptr<TextureMapping2D> mapping, const std::string& filename,
		bool do_trilinear, float max_aniso, ImageWrap wrap_mode, float scale, bool gamma) :
		mapping(std::move(mapping)) {
	mipmap = getTexture(filename, do_trilinear, max_aniso, wrap_mode, scale, gamma);
}

template <typename Tmemory, typename Treturn>
MIPMap<Tmemory>* ImageTexture<Tmemory, Treturn>::getTexture(const std::string& filename, bool do_trilinear, float max_aniso, ImageWrap wrap_mode, float scale, bool gamma) {
	TexInfo info{ filename, do_trilinear, max_aniso, wrap_mode, scale, gamma };
	if (textures.find(info) != textures.end()) {
		return textures[info].get();
	}

	Point2i resolution;
	std::unique_ptr<RGBSpectrum[]> texels = readImage(filename, &resolution);
	MIPMap<Tmemory>* tmp_mipmap = nullptr;
	if (texels) {
		std::unique_ptr<Tmemory[]> converted_texels(new Tmemory[resolution.x * resolution.y]);
		for (int i = 0; i < resolution.x * resolution.y; i++) {
			convertIn(texels[i], &converted_texels[i], scale, gamma);
		}
		tmp_mipmap = new MIPMap<Tmemory>(resolution, converted_texels.get(), do_trilinear, max_aniso, wrap_mode);
	} else {
		Tmemory one = scale;
		tmp_mipmap = new MIPMap<Tmemory>(Point2i(1, 1), &one);
	}
	textures[info].reset(tmp_mipmap);

	return tmp_mipmap;
}

template <typename Tmemory, typename Treturn>
std::map<TexInfo, std::unique_ptr<MIPMap<Tmemory>>>
		ImageTexture<Tmemory, Treturn>::textures;

ImageTexture<float, float>* createImageFloatTexture(const Transform& tex2world, const std::string& filename) {
	// Initialize default parameters
	std::unique_ptr<TextureMapping2D> map;
	map.reset(new UVMapping2D(1., 1., 0., 0.));
	float max_aniso = 8.f;
	bool trilerp = false;
	ImageWrap wrap_mode = ImageWrap::Repeat;
	float scale = 1.f;
	bool gamma = true;
	return new ImageTexture<float, float>(std::move(map), filename, trilerp, max_aniso, wrap_mode, scale, gamma);
}

ImageTexture<RGBSpectrum, Spectrum>* createImageSpectrumTexture(const Transform& tex2world, const std::string& filename) {
	// Initialize default parameters
	std::unique_ptr<TextureMapping2D> map;
	map.reset(new UVMapping2D(1., 1., 0., 0.));
	float max_aniso = 8.f;
	bool trilerp = false;
	ImageWrap wrap_mode = ImageWrap::Repeat;
	float scale = 1.f;
	bool gamma = true;
	return new ImageTexture<RGBSpectrum, Spectrum>(std::move(map), filename, trilerp, max_aniso, wrap_mode, scale, gamma);
}

template class ImageTexture<float, float>;
template class ImageTexture<RGBSpectrum, float>;
template class ImageTexture<RGBSpectrum, Spectrum>;

} // namespace civet