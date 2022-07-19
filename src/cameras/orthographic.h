#ifndef CIVET_ORTHOGRAPHIC_H
#define CIVET_ORTHOGRAPHIC_H

#include <core/camera.h>
#include <core/civet.h>

namespace civet {

class OrthographicCamera : public ProjectiveCamera {
public:
	OrthographicCamera(const AnimatedTransform& ctw, const Bounds2f& screen_window, float open, float close,
			float lens_radius, float focal_distance, Film* f, const Medium* m);

	float generateRay(const CameraSample& sample, Ray* ray) const override;
	float generateRayDifferential(const CameraSample& sample, RayDifferential* rd) const override;

private:
	Vector3f dx_camera, dy_camera;
};

} // namespace civet

#endif // CIVET_ORTHOGRAPHIC_H
