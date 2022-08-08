#ifndef CIVET_SOLID_SHADER_H
#define CIVET_SOLID_SHADER_H

#include <core/shader.h>

namespace civet {

class SolidShader : public Shader {
public:
	SolidShader() :
			Shader("../civet/src/shaders/solid_smooth_vert.glsl", "../civet/src/shaders/solid_smooth_frag.glsl") {}
};

} // namespace civet

#endif // CIVET_SOLID_SHADER_H
