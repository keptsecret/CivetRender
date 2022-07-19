#ifndef CIVET_CAMERA_H
#define CIVET_CAMERA_H

#include <core/civet.h>
#include <core/geometry/transform.h>

namespace civet {

// TODO: remove when sampling is properly implemented
Point2f concentricSamplingDisk(const Point2f& p2f);

struct CameraSample {
	Point2f p_film;
	Point2f p_lens;
	float time;
};

class Camera {
public:
	Camera(const AnimatedTransform& ctw, float open, float close, Film* f, const Medium* m);
	virtual ~Camera();

	virtual float generateRay(const CameraSample& sample, Ray* ray) const = 0;
	virtual float generateRayDifferential(const CameraSample& sample, RayDifferential* rd) const;

	AnimatedTransform camera_to_world;
	const float shutter_open, shutter_close;
	Film* film;
	const Medium* medium;
};

class ProjectiveCamera : public Camera {
public:
	ProjectiveCamera(const AnimatedTransform& ctw, const Transform& cts, const Bounds2f& screen_window, float open, float close,
			float lensr, float focald, Film* f, const Medium* m);

protected:
	Transform camera_to_screen, raster_to_camera;
	Transform screen_to_raster, raster_to_screen;
	float lens_radius, focal_distance;
};

} // namespace civet

#endif // CIVET_CAMERA_H
