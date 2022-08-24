#include <shaders/simple_forward.h>

namespace civet {

void SimpleForwardShader::setPointLights(std::vector<GLPointLight>& lights) {
	for (int i = 0; i < lights.size(); i++) {
		auto diffuse = lights[i].color * lights[i].intensity;

		setBool(("pointLights[" + std::to_string(i) + "].valid").c_str(), lights[i].active && lights[i].cast_shadow);
		setVec3(("pointLights[" + std::to_string(i) + "].position").c_str(), Vector3f(lights[i].position));
		setVec3(("pointLights[" + std::to_string(i) + "].ambient").c_str(), diffuse * 0.1f);
		setVec3(("pointLights[" + std::to_string(i) + "].diffuse").c_str(), diffuse);
		setVec3(("pointLights[" + std::to_string(i) + "].specular").c_str(), 1.0f, 1.0f, 1.0f);
		setFloat(("pointLights[" + std::to_string(i) + "].constant").c_str(), 1.0f);
		setFloat(("pointLights[" + std::to_string(i) + "].linear").c_str(), 0.09f);
		setFloat(("pointLights[" + std::to_string(i) + "].quadratic").c_str(), 0.032f);
	}
}

void SimpleForwardShader::setDirectionalLights(std::vector<GLDirectionalLight>& lights) {
	for (int i = 0; i < 4; i++) {
		auto diffuse = lights[i].color * lights[i].intensity;

		setBool(("dirLights[" + std::to_string(i) + "].valid").c_str(), lights[i].active && lights[i].cast_shadow);
		setVec3(("dirLights[" + std::to_string(i) + "].direction").c_str(), lights[i].direction);
		setVec3(("dirLights[" + std::to_string(i) + "].ambient").c_str(), diffuse * 0.1f);
		setVec3(("dirLights[" + std::to_string(i) + "].diffuse").c_str(), diffuse);
		setVec3(("dirLights[" + std::to_string(i) + "].specular").c_str(), 1.0f, 1.0f, 1.0f);
	}
}

} // namespace civet