#ifndef CIVET_SIMPLE_FORWARD_SHADER_H
#define CIVET_SIMPLE_FORWARD_SHADER_H

#include <core/shader.h>
#include <core/light.h>

namespace civet {

class SimpleForwardShader : public Shader {
public:
	SimpleForwardShader() :
			Shader("../civet/src/shaders/simple_forward_vert.glsl", "../civet/src/shaders/simple_forward_frag.glsl") {}

	void setPointLights(std::vector<GLPointLight>& lights);
	void setDirectionalLights(std::vector<GLDirectionalLight>& lights);
};

} // namespace civet

#endif // CIVET_SIMPLE_FORWARD_SHADER_H
