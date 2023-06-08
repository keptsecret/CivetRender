#ifndef CIVET_SKYBOX_H
#define CIVET_SKYBOX_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/node.h>
#include <core/shader.h>

struct ArHosekSkyModelState;

namespace civet {

struct Atmosphere {
	ArHosekSkyModelState* stateR = nullptr;
	ArHosekSkyModelState* stateG = nullptr;
	ArHosekSkyModelState* stateB = nullptr;

	Vector3f sun_direction;
	float turbidity = 0.f;
	Vector3f albedo;
	float elevation = 0.f;

	~Atmosphere() {
		clear();
	}

	void init(Vector3f sun_dir, Vector3f ground_albedo, float turbidity);

	void clear();
};

struct SkyboxParameters {
	Vector3f sun_direction{ 0, 1, 0 };
	unsigned int resolution = 1024;
	Vector3f ground_color{ 0.5f, 0.5f, 0.5f };
	float turbidity = 2.f;
	float exposure = 16.f;
};

class Skybox : public Node {
public:
	Skybox() :
			Node("Skybox", NodeType::SkyBox, false) {}

	void init(const SkyboxParameters& params);
	void update(const SkyboxParameters& params);
	void draw(const Transform& projection, const Transform& view);
	void renderSkyboxToTexture(const Vector3f& sun_dir);
	Spectrum sampleSky(const Atmosphere& atmosphere, const Vector3f& sample_dir);

	void resetEditingParameters() { editing_params = parameters; }

	SkyboxParameters editing_params;
	Atmosphere atmosphere;

private:
	Vector3f mapToDirection(float x, float y, int s);

	static Vector3f clampSkyExposure(const Vector3f& color, float exposure) {
		return std::exp2(-exposure) * color;
	}

	SkyboxParameters parameters;

	std::vector<std::vector<Vector3f>> skybox_data;
	unsigned int cubemap;

	Shader shader;
	unsigned int VAO;
	unsigned int VBO;
	static const float vertices[3 * 36];
};

} // namespace civet

#endif // CIVET_SKYBOX_H
