#include <cameras/environment.h>

namespace civet {

EnvironmentCamera::EnvironmentCamera(const AnimatedTransform& ctw, float open, float close, Film* f, const Medium* m) :
		Camera(ctw, open, close, f, m) {
}
float EnvironmentCamera::generateRay(const CameraSample& sample, Ray* ray) const {
	float theta = Pi * sample.p_film.y / film->full_resolution.y;
	float phi = 2 * Pi * sample.p_film.x / film->full_resolution.x;
	Vector3f dir(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
	*ray = Ray(Point3f(0, 0, 0), dir, Infinity, lerp(sample.time, shutter_open, shutter_close));
	ray->medium = medium;
	*ray = camera_to_world(*ray);
	return 1;
}

} // namespace civet