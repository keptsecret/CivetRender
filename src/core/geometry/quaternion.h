#ifndef CIVET_QUATERNION_H
#define CIVET_QUATERNION_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>

namespace civet {

class Quaternion {
public:
	Quaternion() :
			v(0, 0, 0), w(1) {}

	Quaternion(const Transform& t);

	Quaternion& operator+=(const Quaternion& q) {
		v += q.v;
		w += q.w;
		return *this;
	}

	friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
		Quaternion ret = q1;
		return ret += q2;
	}

	Quaternion& operator-=(const Quaternion& q) {
		v -= q.v;
		w -= q.w;
		return *this;
	}

	Quaternion operator-() const {
		Quaternion ret;
		ret.v = -v;
		ret.w = -w;
		return ret;
	}

	friend Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
		Quaternion ret = q1;
		return ret -= q2;
	}

	Quaternion& operator*=(float f) {
		v *= f;
		w *= f;
		return *this;
	}

	Quaternion operator*(float f) const {
		Quaternion ret = *this;
		ret.v *= f;
		ret.w *= f;
		return ret;
	}

	Quaternion& operator/=(float f) {
		v /= f;
		w /= f;
		return *this;
	}

	Quaternion operator/(float f) const {
		Quaternion ret = *this;
		ret.v /= f;
		ret.w /= f;
		return ret;
	}

	Transform toTransform() const;

	Vector3f v;
	float w;
};

Quaternion slerp(float t, const Quaternion& q1, const Quaternion& q2);

// inline functions
inline Quaternion operator*(float f, const Quaternion& q) { return q * f; }

inline float dot(const Quaternion& q1, const Quaternion& q2) {
	return dot(q1.v, q2.v) + q1.w * q2.w;
}

inline Quaternion normalize(const Quaternion& q) {
	return q / std::sqrt(dot(q, q));
}

} // namespace civet

#endif // CIVET_QUATERNION_H
