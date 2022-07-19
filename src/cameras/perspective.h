#ifndef CIVET_PERSPECTIVE_H
#define CIVET_PERSPECTIVE_H

#include <core/camera.h>
#include <core/civet.h>

namespace civet {

class PerspectiveCamera : public ProjectiveCamera {
public:
	PerspectiveCamera(const AnimatedTransform& ctw, const Bounds2f& screen_window, float open, float close,
			float lens_radius, float focal_distance, float fov, Film* f, const Medium* m);

	float generateRay(const CameraSample& sample, Ray* ray) const override;
	float generateRayDifferential(const CameraSample& sample, RayDifferential* rd) const override;

private:
	Vector3f dx_camera, dy_camera;
};

} // namespace civet

#endif // CIVET_PERSPECTIVE_H
