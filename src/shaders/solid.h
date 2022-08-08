#ifndef CIVET_SOLID_SHADER_H
#define CIVET_SOLID_SHADER_H

#include <core/shader.h>

namespace civet {

class SolidShader : public Shader {
public:
	SolidShader() :
			Shader("../civet/src/shaders/solid_shader.vert", "../civet/src/shaders/solid_shader.frag") {}
};

} // namespace civet

#endif // CIVET_SOLID_SHADER_H
