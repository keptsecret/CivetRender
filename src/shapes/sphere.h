#ifndef CIVET_SPHERE_H
#define CIVET_SPHERE_H

#include <core/shape.h>

namespace civet {

class Sphere : public Shape {
public:
	CIVET_CPU_GPU
	Sphere(const Transform* otw, const Transform* wto, bool _reverse_orientation, float r, float _z_min, float _z_max, float phi) :
			Shape(otw, wto, _reverse_orientation),
			radius(r),
			z_min(clamp(std::min(_z_min, _z_max), -r, r)),
			z_max(clamp(std::min(_z_min, _z_max), -r, r)),
			theta_min(std::acos(clamp(_z_min / r, -1, 1))),
			theta_max(std::acos(clamp(_z_max / r, -1, 1))),
			phi_max(radians(clamp(phi, 0, 360))) {}

	CIVET_CPU_GPU
	Bounds3f objectBound() const override;

	CIVET_CPU_GPU
	bool intersect(const Ray &ray, float *t_hit, SurfaceInteraction *isect, bool test_alpha_texture = true) const override;

	CIVET_CPU_GPU
	bool intersectP(const Ray &ray, bool test_alpha_texture = true) const override;

	CIVET_CPU_GPU
	float area() const override;

private:
	const float radius;
	const float z_min, z_max;
	const float theta_min, theta_max, phi_max;
};

} // namespace civet

#endif // CIVET_SPHERE_H
