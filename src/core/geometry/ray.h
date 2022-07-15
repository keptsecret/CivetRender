#ifndef CIVET_RAY_H
#define CIVET_RAY_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>

namespace civet {

class Ray {
public:
	CIVET_CPU_GPU
	Ray() :
			t_max(Infinity), time(0.f), medium(nullptr) {}

	CIVET_CPU_GPU
	Ray(const Point3f& origin, const Vector3f& direction, float tm = Infinity, float t = 0.f, const Medium* m = nullptr) :
			o(origin), d(direction), t_max(tm), time(t), medium(m) {
		///< Moved ray-triangle intersect coefficients to be stored in the ray class instead
		// Permute ray so largest ray dimension in z-axis (x and y axis arbitrary) --> +z direction
		isect_coeff.kz = maxDimension(abs(d));
		isect_coeff.kx = isect_coeff.kz + 1;
		if (isect_coeff.kx == 3) {
			isect_coeff.kx = 0;
		}
		isect_coeff.ky = isect_coeff.kx + 1;
		if (isect_coeff.ky == 3) {
			isect_coeff.ky = 0;
		}

		// Shear to align ray along z-axis
		isect_coeff.Sx = -d.x / d.z;
		isect_coeff.Sy = -d.y / d.z;
		isect_coeff.Sz = 1 / d.z;
	}

	CIVET_CPU_GPU
	Point3f operator()(const float t) const { return o + t * d; }

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return o.hasNaNs() || d.hasNaNs() || civet::isNaN(t_max);
	}

	Point3f o;
	Vector3f d;
	mutable float t_max;
	float time;
	const Medium* medium;

	struct {
		int kx = 1, ky = 1, kz = 1;
		float Sx = 0, Sy = 0, Sz = 0;
	} isect_coeff;
};

class RayDifferential : public Ray {
public:
	RayDifferential() :
			has_differentials(false) {}

	CIVET_CPU_GPU
	RayDifferential(const Point3f& origin, const Vector3f& direction, float tm = Infinity, float t = 0.f, const Medium* m = nullptr) :
			Ray(origin, direction, tm, t, m), has_differentials(false) {}

	CIVET_CPU_GPU
	RayDifferential(const Ray& ray) :
			Ray(ray), has_differentials(false) {}

	CIVET_CPU_GPU
	void scaleDifferentials(float s) {
		rx_origin = o + (rx_origin - o) * s;
		ry_origin = o + (ry_origin - o) * s;
		rx_direction = d + (rx_direction - d) * s;
		ry_direction = d + (ry_direction - d) * s;
	}

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return Ray::hasNaNs() || (has_differentials && (rx_origin.hasNaNs() || ry_origin.hasNaNs() || rx_direction.hasNaNs() || ry_direction.hasNaNs()));
	}

	bool has_differentials;
	Point3f rx_origin, ry_origin;
	Vector3f rx_direction, ry_direction;
};

} // namespace civet

#endif // CIVET_RAY_H
