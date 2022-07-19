#ifndef CIVET_ENVIRONMENT_H
#define CIVET_ENVIRONMENT_H

#include <core/civet.h>
#include <core/camera.h>

namespace civet {

class EnvironmentCamera : public Camera {
public:
	EnvironmentCamera(const AnimatedTransform& ctw, float open, float close, Film* f, const Medium* m);

	float generateRay(const CameraSample& sample, Ray* ray) const override;
};

} // namespace civet

#endif // CIVET_ENVIRONMENT_H
