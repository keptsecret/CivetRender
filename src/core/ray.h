#ifndef CIVET_RAY_H
#define CIVET_RAY_H

#include "civet.h"
#include "vecmath.h"

namespace civet {

class Ray {
public:
	Ray() :
			t_max(Infinity), time(0.f), medium(nullptr) {}

	CIVET_CPU_GPU
	Ray(const Point3f& origin, const Vector3f& direction, float tm = Infinity, float t = 0.f, const Medium* m = nullptr) :
			o(origin), d(direction), t_max(tm), time(t), medium(m) {}

	CIVET_CPU_GPU
	Point3f operator()(const float t) const { return o + t * d; }

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return o.hasNaNs() || d.hasNaNs() || std::isnan(t_max);
	}

	friend std::ostream &operator<<(std::ostream &os, const Ray &r) {
		os << "[o=" << r.o << ", d=" << r.d << ", tMax=" << r.t_max
		   << ", time=" << r.time << "]";
		return os;
	}

	Point3f o;
	Vector3f d;
	mutable float t_max;
	float time;
	const Medium* medium;
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

	friend std::ostream &operator<<(std::ostream &os, const RayDifferential &r) {
		os << "[ " << (Ray &)r << " has differentials: " <<
				(r.has_differentials ? "true" : "false") << ", xo = " << r.rx_origin <<
				", xd = " << r.rx_direction << ", yo = " << r.ry_origin << ", yd = " <<
				r.ry_direction;
		return os;
	}

	bool has_differentials;
	Point3f rx_origin, ry_origin;
	Vector3f rx_direction, ry_direction;
};

}


#endif // CIVET_RAY_H
