/*
 * Rayleigh and Mie scattering method for drawing the atmosphere and sky
 * as described by Nishita's paper "Display of the Earth Taking into Account Atmospheric Scattering"
 * Adapted from: https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/simulating-sky/simulating-colors-of-the-sky.html
 */

#ifndef CIVET_SKYBOX_H
#define CIVET_SKYBOX_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/node.h>
#include <core/shader.h>

namespace civet {

struct Atmosphere {
	Vector3f sun_direction{ 0, 1, 0 }; // normalized sun direction
	float earth_radius = 6360e3; // radius to ground (Earth), Re or Rg
	float atmosphere_radius = 6420e3; // atmosphere radius, R or Ra

	float Hr = 7994; // Rayleigh, thickness of atmosphere for uniform density
	float Hm = 1200; // Mie version

	bool computeIncidentLight(const Vector3f& o, const Vector3f& d, float tmin, float tmax,
			Vector3f* out, unsigned int num_samples, unsigned int num_light_samples) const;

	static const Vector3f beta_R;
	static const Vector3f beta_M;
};

struct SkyboxParameters {
	Vector3f sun_direction{ 0, 1, 0 };
	unsigned int resolution = 1024;
	unsigned int samples_per_pixel = 4;
	Vector3f ground_color{ 0.5f, 0.5f, 0.5f };
};

class Skybox : public Node {
public:
	Skybox() :
			Node("Skybox", NodeType::SkyBox, false) {}

	void init(const SkyboxParameters& params);
	void update(const SkyboxParameters& params);
	void draw(const Transform& projection, const Transform& view);
	void renderSkyboxToTexture(const Vector3f& sun_dir);

	void resetEditingParameters() { editing_params = parameters; }

	SkyboxParameters editing_params;

	static const int sun_intensity = 20; // magic number

private:
	void renderSky(int face, float fov, float aspect, int spp);
	Vector3f mapToDirection(float x, float y, int s);

	Atmosphere atmosphere;
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
