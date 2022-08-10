#ifndef CIVET_TRANSFORM_H
#define CIVET_TRANSFORM_H

#include <core/civet.h>
#include <core/geometry/quaternion.h>
#include <core/geometry/ray.h>
#include <core/geometry/vecmath.h>

namespace civet {

struct Matrix4 {
	CIVET_CPU_GPU
	Matrix4() {
		m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
		m[0][1] = m[0][2] = m[0][3] = m[1][0] = m[1][2] = m[1][3] = m[2][0] =
				m[2][1] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0.0f;
	}

	CIVET_CPU_GPU
	Matrix4(float mat[4][4]);

	CIVET_CPU_GPU
	Matrix4(float t00, float t01, float t02, float t03,
			float t10, float t11, float t12, float t13,
			float t20, float t21, float t22, float t23,
			float t30, float t31, float t32, float t33);

	CIVET_CPU_GPU
	bool operator==(const Matrix4& m2) const {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				if (m[i][j] != m2.m[i][j]) {
					return false;
				}
			}
		}
		return true;
	}

	CIVET_CPU_GPU
	bool operator!=(const Matrix4& m2) const {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				if (m[i][j] != m2.m[i][j]) {
					return true;
				}
			}
		}
		return false;
	}

	CIVET_CPU_GPU
	friend Matrix4 transpose(const Matrix4& mat);

	void print(FILE* f) const {
		fprintf(f, "[ ");
		for (int i = 0; i < 4; ++i) {
			fprintf(f, "  [ ");
			for (int j = 0; j < 4; ++j) {
				fprintf(f, "%f", m[i][j]);
				if (j != 3) {
					fprintf(f, ", ");
				}
			}
			fprintf(f, " ]\n");
		}
		fprintf(f, " ] ");
	}

	CIVET_CPU_GPU
	static Matrix4 matmul(const Matrix4& m1, const Matrix4& m2) {
		Matrix4 r;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
						m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
			}
		}
		return r;
	}

	CIVET_CPU_GPU
	friend Matrix4 inverse(const Matrix4& mat);

	friend std::ostream& operator<<(std::ostream& os, const Matrix4& mat) {
		// clang-format off
        os << "[ [ " << mat.m[0][0] << ", " << mat.m[0][1] << ", " << mat.m[0][2] << ", " << mat.m[0][3] << " ] \n"
			<< "  [ " << mat.m[1][0] << ", " << mat.m[1][1] << ", " << mat.m[1][2] << ", " << mat.m[1][3] << " ] \n"
			<< "  [ " << mat.m[2][0] << ", " << mat.m[2][1] << ", " << mat.m[2][2] << ", " << mat.m[2][3] << " ] \n"
			<< "  [ " << mat.m[3][0] << ", " << mat.m[3][1] << ", " << mat.m[3][2] << ", " << mat.m[3][3] << " ] ]";

		// clang-format on
		return os;
	}

	float m[4][4];
};

class Transform {
public:
	Transform() {}

	CIVET_CPU_GPU
	Transform(const float mat[4][4]) {
		m = Matrix4(mat[0][0], mat[0][1], mat[0][2], mat[0][3],
				mat[1][0], mat[1][1], mat[1][2], mat[1][3],
				mat[2][0], mat[2][1], mat[2][2], mat[2][3],
				mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
		m_inv = inverse(m);
	}

	CIVET_CPU_GPU
	Transform(const Matrix4& mat) :
			m(mat), m_inv(inverse(mat)) {}

	CIVET_CPU_GPU
	Transform(const Matrix4& mat, const Matrix4& mat_inv) :
			m(mat), m_inv(mat_inv) {}

	CIVET_CPU_GPU
	friend Transform inverse(const Transform& t) {
		return Transform(t.m_inv, t.m);
	}

	CIVET_CPU_GPU
	friend Transform transpose(const Transform& t) {
		return Transform(transpose(t.m), transpose(t.m_inv));
	}

	CIVET_CPU_GPU
	bool operator==(const Transform& t) const {
		return t.m == m && t.m_inv == m_inv;
	}

	CIVET_CPU_GPU
	bool operator!=(const Transform& t) const {
		return t.m != m || t.m_inv != m_inv;
	}

	CIVET_CPU_GPU
	bool operator<(const Transform& t2) const {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				if (m.m[i][j] < t2.m.m[i][j]) {
					return true;
				}
				if (m.m[i][j] > t2.m.m[i][j]) {
					return false;
				}
			}
		}
		return false;
	}

	CIVET_CPU_GPU
	bool isIdentity() const {
		return (m.m[0][0] == 1.f && m.m[0][1] == 0.f && m.m[0][2] == 0.f && m.m[0][3] == 0.f &&
				m.m[1][0] == 0.f && m.m[1][1] == 1.f && m.m[1][2] == 0.f && m.m[1][3] == 0.f &&
				m.m[2][0] == 0.f && m.m[2][1] == 0.f && m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
				m.m[3][0] == 0.f && m.m[3][1] == 0.f && m.m[3][2] == 0.f && m.m[3][3] == 1.f);
	}

	const Matrix4& getMatrix() const { return m; }
	const Matrix4& getInverseMatrix() const { return m_inv; }

	CIVET_CPU_GPU
	bool hasScale() const {
		float la2 = (*this)(Vector3f(1, 0, 0)).lengthSquared();
		float lb2 = (*this)(Vector3f(0, 1, 0)).lengthSquared();
		float lc2 = (*this)(Vector3f(0, 0, 1)).lengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
		return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
	}

	// Apply transforms
	template <typename T>
	CIVET_CPU_GPU inline Point3<T> operator()(const Point3<T>& p) const;
	template <typename T>
	CIVET_CPU_GPU inline Vector3<T> operator()(const Vector3<T>& v) const;
	template <typename T>
	CIVET_CPU_GPU inline Normal3<T> operator()(const Normal3<T>&) const;
	CIVET_CPU_GPU
	inline Ray operator()(const Ray& r) const;
	CIVET_CPU_GPU
	inline RayDifferential operator()(const RayDifferential& r) const;

	CIVET_CPU_GPU
	Bounds3f operator()(const Bounds3f& b) const;
	CIVET_CPU_GPU
	Transform operator*(const Transform& t2) const;

	CIVET_CPU_GPU
	bool swapsHandedness() const;

	CIVET_CPU_GPU
	SurfaceInteraction operator()(const SurfaceInteraction& si) const;

	// Apply transforms with rounding errors
	// TODO: implement rounding error transforms
	template <typename T>
	CIVET_CPU_GPU inline Point3<T> operator()(const Point3<T>& pt, Vector3<T>* p_error) const;
	template <typename T>
	CIVET_CPU_GPU inline Point3<T> operator()(const Point3<T>& pt, const Vector3<T>& pt_error, Vector3<T>* abs_error) const;
	//	template <typename T>
	//	CIVET_CPU_GPU inline Vector3<T> operator()(const Vector3<T>& v, Vector3<T>* vTransError) const;
	//	template <typename T>
	//	CIVET_CPU_GPU inline Vector3<T> operator()(const Vector3<T>& v, const Vector3<T>& vError, Vector3<T>* vTransError) const;
	//	CIVET_CPU_GPU inline Ray operator()(const Ray& r, Vector3f* oError, Vector3f* dError) const;
	//	CIVET_CPU_GPU inline Ray operator()(const Ray& r, const Vector3f& oErrorIn, const Vector3f& dErrorIn, Vector3f* oErrorOut, Vector3f* dErrorOut) const;

	Matrix4 m, m_inv;
};

CIVET_CPU_GPU
Transform translate(const Vector3f& delta);
CIVET_CPU_GPU
Transform scale(const float x, const float y, const float z);
CIVET_CPU_GPU
Transform rotateX(float theta);
CIVET_CPU_GPU
Transform rotateY(float theta);
CIVET_CPU_GPU
Transform rotateZ(float theta);
CIVET_CPU_GPU
Transform rotate(float theta, Vector3f& axis);

/**
 * Generate a look-at matrix in left-handed coordinates\n
 * For use with the PATH-TRACER
 * @param position camera position
 * @param target target the camera is facing
 * @param up local up vector of camera, usually (0, 1, 0)
 * @return a Transform view matrix in left-handed coordinates
 */
CIVET_CPU_GPU
Transform lookAtLH(const Point3f& position, const Point3f& target, const Vector3f& up);
/**
* Generate a look-at matrix in right-handed coordinates\n
* Used primarily with OpenGL coordinate system
* @param position camera position
* @param target target the camera is facing
* @param up local up vector of camera, usually (0, 1, 0)
* @return a Transform view matrix in right-handed coordinates
 */
CIVET_CPU_GPU
Transform lookAtRH(const Point3f& position, const Point3f& target, const Vector3f& up);

CIVET_CPU_GPU
Transform orthographic(float z_near, float z_far);
CIVET_CPU_GPU
Transform orthographic(float left, float right, float bottom, float top, float z_near, float z_far);
CIVET_CPU_GPU
Transform perspective(float fov, float n, float f);
CIVET_CPU_GPU
Transform perspective(float fov, float aspect, float n, float f);

// Transform Inline Functions
template <typename T>
CIVET_CPU_GPU inline Point3<T> Transform::operator()(const Point3<T>& p) const {
	T x = p.x, y = p.y, z = p.z;
	T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
	T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
	T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
	T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
	if (wp == 1) {
		return Point3<T>(xp, yp, zp);
	} else {
		return Point3<T>(xp, yp, zp) / wp;
	}
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> Transform::operator()(const Vector3<T>& v) const {
	T x = v.x, y = v.y, z = v.z;
	return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
			m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
			m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> Transform::operator()(const Normal3<T>& n) const {
	T x = n.x, y = n.y, z = n.z;
	return Normal3<T>(m_inv.m[0][0] * x + m_inv.m[1][0] * y + m_inv.m[2][0] * z,
			m_inv.m[0][1] * x + m_inv.m[1][1] * y + m_inv.m[2][1] * z,
			m_inv.m[0][2] * x + m_inv.m[1][2] * y + m_inv.m[2][2] * z);
}

CIVET_CPU_GPU
inline Ray Transform::operator()(const Ray& r) const {
	Vector3f oError;
	Point3f o = (*this)(r.o, &oError);
	Vector3f d = (*this)(r.d);

	// Offset ray origin to edge of error bounds and compute t_max
	float lengthSquared = d.lengthSquared();
	float tMax = r.t_max;
	if (lengthSquared > 0) {
		float dt = dot(abs(d), oError) / lengthSquared;
		o += d * dt;
		tMax -= dt;
	}
	return Ray(o, d, tMax, r.time, r.medium);
}

CIVET_CPU_GPU
inline RayDifferential Transform::operator()(const RayDifferential& r) const {
	Ray tr = (*this)(Ray(r));
	RayDifferential ret(tr.o, tr.d, tr.t_max, tr.time, tr.medium);
	ret.has_differentials = r.has_differentials;
	ret.ry_origin = (*this)(r.ry_origin);
	ret.ry_origin = (*this)(r.ry_origin);
	ret.ry_direction = (*this)(r.ry_direction);
	ret.ry_direction = (*this)(r.ry_direction);
	return ret;
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> Transform::operator()(const Point3<T>& p, Vector3<T>* p_error) const {
	T x = p.x, y = p.y, z = p.z;
	// Compute transformed coordinates from point _pt_
	T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
	T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
	T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
	T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

	// Compute absolute error for transformed point
	T xAbsSum = (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
			std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
	T yAbsSum = (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
			std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
	T zAbsSum = (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
			std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
	*p_error = gamma(3) * Vector3<T>(xAbsSum, yAbsSum, zAbsSum);
	if (wp == 1) {
		return Point3<T>(xp, yp, zp);
	} else {
		return Point3<T>(xp, yp, zp) / wp;
	}
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> Transform::operator()(const Point3<T>& pt, const Vector3<T>& pt_error, Vector3<T>* abs_error) const {
	T x = pt.x, y = pt.y, z = pt.z;
	T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
	T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
	T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
	T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
	abs_error->x =
			(gamma(3) + (T)1) *
					(std::abs(m.m[0][0]) * pt_error.x + std::abs(m.m[0][1]) * pt_error.y +
							std::abs(m.m[0][2]) * pt_error.z) +
			gamma(3) * (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) + std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
	abs_error->y =
			(gamma(3) + (T)1) *
					(std::abs(m.m[1][0]) * pt_error.x + std::abs(m.m[1][1]) * pt_error.y +
							std::abs(m.m[1][2]) * pt_error.z) +
			gamma(3) * (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) + std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
	abs_error->z =
			(gamma(3) + (T)1) *
					(std::abs(m.m[2][0]) * pt_error.x + std::abs(m.m[2][1]) * pt_error.y +
							std::abs(m.m[2][2]) * pt_error.z) +
			gamma(3) * (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) + std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
	if (wp == 1.) {
		return Point3<T>(xp, yp, zp);
	} else {
		return Point3<T>(xp, yp, zp) / wp;
	}
}

class AnimatedTransform {
public:
	AnimatedTransform(const Transform* start_tr, float start_t, const Transform* end_tr, float end_t);

	void decompose(const Matrix4& m, Vector3f* T, Quaternion* Rquat, Matrix4* S);

	void interpolate(float time, Transform* t) const;

	Ray operator()(const Ray& r) const;
	RayDifferential operator()(const RayDifferential& r) const;
	Point3f operator()(float time, const Point3f& p) const;
	Vector3f operator()(float time, const Vector3f& v) const;

	Bounds3f motionBounds(const Bounds3f& b) const;
	Bounds3f boundPointMotion(const Point3f& p) const;

private:
	const Transform *start_transform, *end_transform;
	const float start_time, end_time;
	const bool is_animated;
	Vector3f T[2];
	Quaternion R[2];
	Matrix4 S[2];
	bool has_rotation;

	struct DerivativeTerm {
		DerivativeTerm() {}
		DerivativeTerm(float c, float x, float y, float z) :
				kc(c), kx(x), ky(y), kz(z) {}

		float kc, kx, ky, kz;
		float eval(const Point3f& p) const {
			return kc + kx * p.x + ky * p.y + kz * p.z;
		}
	};
	DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
};

} // namespace civet

#endif // CIVET_TRANSFORM_H
