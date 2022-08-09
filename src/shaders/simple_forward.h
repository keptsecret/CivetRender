#ifndef CIVET_SIMPLE_FORWARD_SHADER_H
#define CIVET_SIMPLE_FORWARD_SHADER_H

#include <core/shader.h>

namespace civet {

class SimpleForwardShader : public Shader {
public:
	SimpleForwardShader() :
			Shader("../civet/src/shaders/simple_forward_vert.glsl", "../civet/src/shaders/simple_forward_frag.glsl") {}

	void use(std::vector<Point3f>& light_positions) override {
		Shader::use(light_positions);

		for (int i = 0; i < light_positions.size(); i++) {
			setBool(("pointLights[" + std::to_string(i) + "].valid").c_str(), true);
			setVec3(("pointLights[" + std::to_string(i) + "].position").c_str(), Vector3f(light_positions[i]));
			setVec3(("pointLights[" + std::to_string(i) + "].ambient").c_str(), 0.05f, 0.05f, 0.05f);
			setVec3(("pointLights[" + std::to_string(i) + "].diffuse").c_str(), 0.8f, 0.8f, 0.8f);
			setVec3(("pointLights[" + std::to_string(i) + "].specular").c_str(), 1.0f, 1.0f, 1.0f);
			setFloat(("pointLights[" + std::to_string(i) + "].constant").c_str(), 1.0f);
			setFloat(("pointLights[" + std::to_string(i) + "].linear").c_str(), 0.09f);
			setFloat(("pointLights[" + std::to_string(i) + "].quadratic").c_str(), 0.032f);
		}

		for (int i = light_positions.size(); i < 32; i++) {
			setBool(("pointLights[" + std::to_string(i) + "].valid").c_str(), false);
		}
	}
};

} // namespace civet

#endif // CIVET_SIMPLE_FORWARD_SHADER_H
