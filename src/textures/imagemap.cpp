#include <textures/imagemap.h>

#include <utils/imageio.h>

namespace civet {

template <typename Tmemory, typename Treturn>
ImageTexture<Tmemory, Treturn>::ImageTexture(std::unique_ptr<TextureMapping2D> mapping, const std::string& filename,
		bool do_trilinear, float max_aniso, ImageWrap wrap_mode, float scale, bool gamma) : mapping(std::move(mapping)) {
	mipmap = getTexture(filename, do_trilinear, max_aniso, wrap_mode, scale, gamma);
}

template <typename Tmemory, typename Treturn>
MIPMap<Tmemory>* ImageTexture<Tmemory, Treturn>::getTexture(const std::string& filename, bool do_trilinear,
		float max_aniso, ImageWrap wrap_mode, float scale, bool gamma) {
	TexInfo info{filename, do_trilinear, max_aniso, wrap_mode, scale, gamma};
	if (textures.find(info) != textures.end()) {
		return textures[info].get();
	}

	Point2i resolution;
	std::unique_ptr<RGBSpectrum[]> texels = readImage(filename, &resolution);
	MIPMap<Tmemory>* tmp_mipmap = nullptr;
	if (texels) {
		std::unique_ptr<Tmemory> converted_texels(new Tmemory[resolution.x * resolution.y]);
		for (int i = 0; i < resolution.x * resolution.y; i++) {
			convertIn(texels[i], &converted_texels[i], scale, gamma);
		}
		tmp_mipmap = new MIPMap<Tmemory>(resolution, converted_texels.get(), do_trilinear, max_aniso, gamma);
	} else {
		Tmemory one = scale;
		tmp_mipmap = new MIPMap<Tmemory>(Point2i(1,1), &one);
	}
	textures[info].reset(tmp_mipmap);

	return tmp_mipmap;
}

} // namespace civet