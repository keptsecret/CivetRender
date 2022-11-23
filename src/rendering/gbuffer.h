#ifndef CIVET_GBUFFER_H
#define CIVET_GBUFFER_H

#include <core/civet.h>

namespace civet {

class GBuffer {
public:
	enum GBUFFER_TEXTURE_TYPE {
		GBUFFER_TEXTURE_POSITION,
		GBUFFER_TEXTURE_DIFFUSE,
		GBUFFER_TEXTURE_SPECULAR,
		GBUFFER_TEXTURE_NORMAL,
		GBUFFER_NUM_TEXTURES
	};

	GBuffer();
	~GBuffer();

	bool init(unsigned int width, unsigned int height);

	void bindWrite();
	void bindRead();

private:
	unsigned int FBO;
	unsigned int textures[GBUFFER_NUM_TEXTURES];
	unsigned int num_textures = GBUFFER_NUM_TEXTURES;
	unsigned int depth_map;
};

} // namespace civet

#endif // CIVET_GBUFFER_H
