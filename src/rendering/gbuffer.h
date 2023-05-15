#ifndef CIVET_GBUFFER_H
#define CIVET_GBUFFER_H

#include <core/civet.h>

namespace civet {

class GBuffer {
public:
	enum GBUFFER_TEXTURE_TYPE {
		GBUFFER_TEXTURE_POSITION,
		GBUFFER_TEXTURE_ALBEDO,
		GBUFFER_TEXTURE_METALLICROUGHAO,
		GBUFFER_TEXTURE_NORMAL,
		GBUFFER_NUM_TEXTURES
	};

	GBuffer();
	~GBuffer();

	bool init(unsigned int width, unsigned int height);
	void start();

	void bindGeomPass();
	void bindStencilPass();
	void bindLightPass();
	void bindPostProcessPass();
	void bindFinalPass();

	unsigned int num_textures = GBUFFER_NUM_TEXTURES;

private:
	unsigned int FBO;
	unsigned int textures[GBUFFER_NUM_TEXTURES];
	unsigned int depth_map;

	unsigned int raw_texture;
	unsigned int final_texture;
};

} // namespace civet

#endif // CIVET_GBUFFER_H
