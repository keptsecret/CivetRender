#include <core/geometry/transform.h>

#include <core/interaction.h>

namespace civet {

bool solveLinearSystem(const float A[2][2], const float B[2], float* x0, float* x1) {
	float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
	if (std::abs(det) < 1e-10f) {
		return false;
	}
	*x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
	*x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
	if (civet::isNaN(*x0) || civet::isNaN(*x1)) {
		return false;
	}
	return true;
}

CIVET_CPU_GPU
Matrix4::Matrix4(float mat[4][4]) {
	memcpy(m, mat, 16 * sizeof(float));
}

CIVET_CPU_GPU
Matrix4::Matrix4(float t00, float t01, float t02, float t03,
		float t10, float t11, float t12, float t13,
		float t20, float t21, float t22, float t23,
		float t30, float t31, float t32, float t33) {
	m[0][0] = t00;
	m[0][1] = t01;
	m[0][2] = t02;
	m[0][3] = t03;
	m[1][0] = t10;
	m[1][1] = t11;
	m[1][2] = t12;
	m[1][3] = t13;
	m[2][0] = t20;
	m[2][1] = t21;
	m[2][2] = t22;
	m[2][3] = t23;
	m[3][0] = t30;
	m[3][1] = t31;
	m[3][2] = t32;
	m[3][3] = t33;
}

CIVET_CPU_GPU
Matrix4 transpose(const Matrix4& mat) {
	return Matrix4(mat.m[0][0], mat.m[1][0], mat.m[2][0], mat.m[3][0],
			mat.m[0][1], mat.m[1][1], mat.m[2][1], mat.m[3][1],
			mat.m[0][2], mat.m[1][2], mat.m[2][2], mat.m[3][2],
			mat.m[0][3], mat.m[1][3], mat.m[2][3], mat.m[3][3]);
}

CIVET_CPU_GPU
Matrix4 inverse(const Matrix4& mat) {
	int indxc[4], indxr[4];
	int ipiv[4] = { 0, 0, 0, 0 };
	float minv[4][4];
	memcpy(minv, mat.m, 4 * 4 * sizeof(float));
	for (int i = 0; i < 4; i++) {
		int irow = 0, icol = 0;
		float big = 0.f;
		// Choose pivot
		for (int j = 0; j < 4; j++) {
			if (ipiv[j] != 1) {
				for (int k = 0; k < 4; k++) {
					if (ipiv[k] == 0) {
						if (std::abs(minv[j][k]) >= big) {
							big = float(std::abs(minv[j][k]));
							irow = j;
							icol = k;
						}
					} else if (ipiv[k] > 1) {
						printf("ERROR::Matrix4: Singular matrix in MatrixInvert\n");
					}
				}
			}
		}
		++ipiv[icol];
		// Swap rows _irow_ and _icol_ for pivot
		if (irow != icol) {
			for (int k = 0; k < 4; ++k) {
				swapElem(minv[irow][k], minv[icol][k]);
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (minv[icol][icol] == 0.f) {
			printf("ERROR::Matrix4: Singular matrix in MatrixInvert\n");
		}

		// Set $mat[icol][icol]$ to one by scaling row _icol_ appropriately
		float pivinv = 1. / minv[icol][icol];
		minv[icol][icol] = 1.;
		for (int j = 0; j < 4; j++) {
			minv[icol][j] *= pivinv;
		}

		// Subtract this row from others to zero out their columns
		for (int j = 0; j < 4; j++) {
			if (j != icol) {
				float save = minv[j][icol];
				minv[j][icol] = 0;
				for (int k = 0; k < 4; k++) {
					minv[j][k] -= minv[icol][k] * save;
				}
			}
		}
	}
	// Swap columns to reflect permutation
	for (int j = 3; j >= 0; j--) {
		if (indxr[j] != indxc[j]) {
			for (int k = 0; k < 4; k++) {
				swapElem(minv[k][indxr[j]], minv[k][indxc[j]]);
			}
		}
	}
	return Matrix4(minv);
}

CIVET_CPU_GPU
Transform translate(const Vector3f& delta) {
	Matrix4 m(1, 0, 0, delta.x,
			0, 1, 0, delta.y,
			0, 0, 1, delta.z,
			0, 0, 0, 1);
	Matrix4 m_inv(1, 0, 0, -delta.x,
			0, 1, 0, -delta.y,
			0, 0, 1, -delta.z,
			0, 0, 0, 1);
	return Transform(m, m_inv);
}

CIVET_CPU_GPU
Transform scale(const float x, const float y, const float z) {
	Matrix4 m(x, 0, 0, 0,
			0, y, 0, 0,
			0, 0, z, 0,
			0, 0, 0, 1);
	Matrix4 m_inv(1 / x, 0, 0, 0,
			0, 1 / y, 0, 0,
			0, 0, 1 / z, 0,
			0, 0, 0, 1);
	return Transform(m, m_inv);
}

CIVET_CPU_GPU
Transform rotateX(float theta) {
	float sin_theta = std::sin(radians(theta));
	float cos_theta = std::cos(radians(theta));
	Matrix4 m(1, 0, 0, 0,
			0, cos_theta, -sin_theta, 0,
			0, sin_theta, cos_theta, 0,
			0, 0, 0, 1);
	return Transform(m, transpose(m));
}

CIVET_CPU_GPU
Transform rotateY(float theta) {
	float sin_theta = std::sin(radians(theta));
	float cos_theta = std::cos(radians(theta));
	Matrix4 m(cos_theta, 0, sin_theta, 0,
			0, 1, 0, 0,
			-sin_theta, 0, cos_theta, 0,
			0, 0, 0, 1);
	return Transform(m, transpose(m));
}

CIVET_CPU_GPU
Transform rotateZ(float theta) {
	float sin_theta = std::sin(radians(theta));
	float cos_theta = std::cos(radians(theta));
	Matrix4 m(cos_theta, -sin_theta, 0, 0,
			0, sin_theta, cos_theta, 0,
			0, 0, 1, 0,
			0, 0, 0, 1);
	return Transform(m, transpose(m));
}

CIVET_CPU_GPU
Transform rotate(float theta, Vector3f& axis) {
	Vector3f a = normalize(axis);
	float sin_theta = std::sin(radians(theta));
	float cos_theta = std::cos(radians(theta));
	Matrix4 m;

	// rotation of first basis vector
	m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cos_theta;
	m.m[0][1] = a.x * a.y + (1 - cos_theta) - a.z * sin_theta;
	m.m[0][2] = a.x * a.z + (1 - cos_theta) - a.y * sin_theta;
	m.m[0][3] = 0;

	// rotation of second basis vector
	m.m[1][0] = a.x * a.y + (1 - cos_theta) - a.x * sin_theta;
	m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cos_theta;
	m.m[1][2] = a.y * a.z + (1 - cos_theta) - a.z * sin_theta;
	m.m[1][3] = 0;

	// rotation of third basis vector
	m.m[2][0] = a.x * a.z + (1 - cos_theta) - a.x * sin_theta;
	m.m[2][1] = a.y * a.z + (1 - cos_theta) - a.y * sin_theta;
	m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cos_theta;
	m.m[2][3] = 0;

	return Transform(m, transpose(m));
}

CIVET_CPU_GPU
Transform lookAtLH(const Point3f& position, const Point3f& target, const Vector3f& up) {
	Matrix4 camera_to_world;
	camera_to_world.m[0][3] = position.x;
	camera_to_world.m[1][3] = position.y;
	camera_to_world.m[2][3] = position.z;
	camera_to_world.m[3][3] = 1;

	Vector3f dir = normalize(target - position);
	Vector3f right = normalize(cross(normalize(up), dir));
	Vector3f new_up = cross(dir, right);

	camera_to_world.m[0][0] = right.x;
	camera_to_world.m[1][0] = right.y;
	camera_to_world.m[2][0] = right.z;
	camera_to_world.m[3][0] = 0;

	camera_to_world.m[0][1] = new_up.x;
	camera_to_world.m[1][1] = new_up.y;
	camera_to_world.m[2][1] = new_up.z;
	camera_to_world.m[3][1] = 0;

	camera_to_world.m[0][2] = dir.x;
	camera_to_world.m[1][2] = dir.y;
	camera_to_world.m[2][2] = dir.z;
	camera_to_world.m[3][2] = 0;

	return Transform(inverse(camera_to_world), camera_to_world);
}

CIVET_CPU_GPU
Transform lookAtRH(const Point3f& position, const Point3f& target, const Vector3f& up) {
	Vector3f dir = normalize(target - position);
	Vector3f right = normalize(cross(dir, normalize(up)));
	Vector3f new_up = cross(right, dir);

	Vector3f eye = Vector3f(position);
	Matrix4 look_at = { right.x, right.y, right.z, -dot(right, eye),
		new_up.x, new_up.y, new_up.z, -dot(new_up, eye),
		-dir.x, -dir.y, -dir.z, dot(dir, eye),
		0, 0, 0, 1 };

	return Transform(look_at, inverse(look_at));
}

CIVET_CPU_GPU
Transform orthographic(float z_near, float z_far) {
	return scale(1, 1, 1 / (z_far - z_near)) * translate(Vector3f(0, 0, -z_near));
}

CIVET_CPU_GPU
Transform orthographic(float left, float right, float bottom, float top, float z_near, float z_far) {
	return scale(2.f / (right - left), 2.f / (top - bottom), 2.f / (z_far - z_near)) *
			translate(Vector3f(-0.5 * (left + right), -0.5 * (top + bottom), -0.5 * (z_near + z_far)));
}

CIVET_CPU_GPU
Transform perspective(float fov, float n, float f) {
	Matrix4 persp(1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, f / (f - n), -f * n / (f - n),
			0, 0, 1, 0);

	float inv_tan_ang = 1 / std::tan(radians(fov) / 2);
	return scale(inv_tan_ang, inv_tan_ang, 1) * Transform(persp);
}

CIVET_CPU_GPU
Transform perspective(float fov, float aspect, float n, float f) {
	float inv_tan_ang = 1 / std::tan(radians(fov) / 2);
	Matrix4 persp(inv_tan_ang / aspect, 0, 0, 0,
			0, inv_tan_ang, 0, 0,
			0, 0, -(f + n) / (f - n), -(2 * f * n) / (f - n),
			0, 0, -1, 0);

	return Transform(persp);
}

CIVET_CPU_GPU
Bounds3f Transform::operator()(const Bounds3f& b) const {
	const Transform& M = *this;
	Bounds3f ret(M(Point3f(b.p_min.x, b.p_min.y, b.p_min.z)));
	ret = bUnion(ret, M(Point3f(b.p_max.x, b.p_min.y, b.p_min.z)));
	ret = bUnion(ret, M(Point3f(b.p_min.x, b.p_max.y, b.p_min.z)));
	ret = bUnion(ret, M(Point3f(b.p_min.x, b.p_min.y, b.p_max.z)));
	ret = bUnion(ret, M(Point3f(b.p_max.x, b.p_max.y, b.p_min.z)));
	ret = bUnion(ret, M(Point3f(b.p_min.x, b.p_max.y, b.p_max.z)));
	ret = bUnion(ret, M(Point3f(b.p_max.x, b.p_min.y, b.p_max.z)));
	ret = bUnion(ret, M(Point3f(b.p_max.x, b.p_max.y, b.p_max.z)));
	return ret;
}

CIVET_CPU_GPU
Transform Transform::operator*(const Transform& t2) const {
	return Transform(Matrix4::matmul(m, t2.m), Matrix4::matmul(t2.m_inv, m_inv));
}

CIVET_CPU_GPU
bool Transform::swapsHandedness() const {
	float det = m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
			m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
			m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
	return det < 0;
}

CIVET_CPU_GPU
SurfaceInteraction Transform::operator()(const SurfaceInteraction& si) const {
	SurfaceInteraction ret;
	const Transform& t = *this;
	ret.p = t(si.p, si.p_error, &ret.p_error);
	ret.n = normalize(t(si.n));
	ret.wo = t(si.wo);
	ret.time = si.time;
	ret.medium_interface = si.medium_interface;
	ret.uv = si.uv;
	ret.shape = si.shape;
	ret.dpdu = t(si.dpdu);
	ret.dpdv = t(si.dpdv);
	ret.dndu = t(si.dndu);
	ret.dndv = t(si.dndv);
	ret.shading.n = normalize(t(si.shading.n));
	ret.shading.dpdu = t(si.shading.dpdu);
	ret.shading.dpdv = t(si.shading.dpdv);
	ret.shading.dndu = t(si.shading.dndu);
	ret.shading.dndv = t(si.shading.dndv);
	ret.dudx = si.dudx;
	ret.dvdx = si.dvdx;
	ret.dudy = si.dudy;
	ret.dvdy = si.dvdy;
	ret.dpdx = t(si.dpdx);
	ret.dpdy = t(si.dpdy);
	ret.bsdf = si.bsdf;
	ret.bssrdf = si.bssrdf;
	ret.primitive = si.primitive;
	// ret.n = faceforward(ret.n, ret.shading.n);
	ret.shading.n = faceforward(ret.shading.n, ret.n);
	return ret;
}

// Class representing intervals of real numbers
class Interval {
public:
	Interval(float v) :
			low(v), high(v) {}
	Interval(float v0, float v1) :
			low(std::min(v0, v1)), high(std::max(v0, v1)) {}

	Interval operator+(const Interval& i) const {
		return Interval(low + i.low, high + i.high);
	}

	Interval operator-(const Interval& i) const {
		return Interval(low - i.high, high - i.low);
	}

	Interval operator*(const Interval& i) const {
		return Interval(std::min(std::min(low * i.low, high * i.low),
								std::min(low * i.high, high * i.high)),
				std::max(std::max(low * i.low, high * i.low),
						std::max(low * i.high, high * i.high)));
	}

	float low, high;
};

inline Interval sin(const Interval& i) {
	float sinLow = std::sin(i.low), sinHigh = std::sin(i.high);
	if (sinLow > sinHigh) {
		swapElem(sinLow, sinHigh);
	}
	if (i.low < Pi / 2 && i.high > Pi / 2) {
		sinHigh = 1.;
	}
	if (i.low < (3.f / 2.f) * Pi && i.high > (3.f / 2.f) * Pi) {
		sinLow = -1.;
	}
	return Interval(sinLow, sinHigh);
}

inline Interval cos(const Interval& i) {
	float cosLow = std::cos(i.low), cosHigh = std::cos(i.high);
	if (cosLow > cosHigh) {
		swapElem(cosLow, cosHigh);
	}
	if (i.low < Pi && i.high > Pi) {
		cosLow = -1.;
	}
	return Interval(cosLow, cosHigh);
}

void intervalFindZeros(float c1, float c2, float c3, float c4, float c5,
		float theta, Interval t_interval, float* zeros, int* zero_count, int depth = 0) {
	Interval range = Interval(c1) +
			(Interval(c2) + Interval(c3) * t_interval) * cos(Interval(2 * theta) * t_interval) +
			(Interval(c4) + Interval(c5) * t_interval) * sin(Interval(2 * theta) * t_interval);
	if (range.low > 0 || range.high < 0 || range.low == range.high) {
		return;
	}

	if (depth > 0) {
		float mid = (t_interval.low + t_interval.high) * 0.5f;
		intervalFindZeros(c1, c2, c3, c4, c5, theta, Interval(t_interval.low, mid), zeros, zero_count, depth - 1);
		intervalFindZeros(c1, c2, c3, c4, c5, theta, Interval(mid, t_interval.high), zeros, zero_count, depth - 1);
	} else {
		float t_newton = (t_interval.low + t_interval.high) * 0.5f;
		for (int i = 0; i < 4; ++i) {
			float f_newton = c1 +
					(c2 + c3 * t_newton) * std::cos(2.f * theta * t_newton) +
					(c4 + c5 * t_newton) * std::sin(2.f * theta * t_newton);
			float f_prime_newton =
					(c3 + 2 * (c4 + c5 * t_newton) * theta) *
							std::cos(2.f * t_newton * theta) +
					(c5 - 2 * (c2 + c3 * t_newton) * theta) *
							std::sin(2.f * t_newton * theta);
			if (f_newton == 0 || f_prime_newton == 0) {
				break;
			}
			t_newton = t_newton - f_newton / f_prime_newton;
		}
		zeros[*zero_count] = t_newton;
		(*zero_count)++;
	}
}

AnimatedTransform::AnimatedTransform(const Transform* start_tr, float start_t, const Transform* end_tr, float end_t) :
		start_transform(start_tr), end_transform(end_tr), start_time(start_t), end_time(end_t), is_animated(*start_tr != *end_tr) {
	decompose(start_tr->m, &T[0], &R[0], &S[0]);
	decompose(end_tr->m, &T[1], &R[1], &S[1]);

	// ensure shortest path is taken
	if (dot(R[0], R[1]) < 0) {
		R[1] = -R[1];
	}

	has_rotation = dot(R[0], R[1]) < 0.9995f;

	// compute terms of motion derivative function
	if (has_rotation) {
		float cos_theta = dot(R[0], R[1]);
		float theta = std::acos(clamp(cos_theta, -1, 1));
		Quaternion qperp = normalize(R[1] - R[0] * cos_theta);

		float t0x = T[0].x;
		float t0y = T[0].y;
		float t0z = T[0].z;
		float t1x = T[1].x;
		float t1y = T[1].y;
		float t1z = T[1].z;
		float q0x = R[0].v.x;
		float q0y = R[0].v.y;
		float q0z = R[0].v.z;
		float q0w = R[0].w;
		float qperpx = qperp.v.x;
		float qperpy = qperp.v.y;
		float qperpz = qperp.v.z;
		float qperpw = qperp.w;
		float s000 = S[0].m[0][0];
		float s001 = S[0].m[0][1];
		float s002 = S[0].m[0][2];
		float s010 = S[0].m[1][0];
		float s011 = S[0].m[1][1];
		float s012 = S[0].m[1][2];
		float s020 = S[0].m[2][0];
		float s021 = S[0].m[2][1];
		float s022 = S[0].m[2][2];
		float s100 = S[1].m[0][0];
		float s101 = S[1].m[0][1];
		float s102 = S[1].m[0][2];
		float s110 = S[1].m[1][0];
		float s111 = S[1].m[1][1];
		float s112 = S[1].m[1][2];
		float s120 = S[1].m[2][0];
		float s121 = S[1].m[2][1];
		float s122 = S[1].m[2][2];

		c1[0] = DerivativeTerm(
				-t0x + t1x,
				(-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
								s000 +
						q0w * q0z * s010 - qperpx * qperpy * s010 +
						qperpw * qperpz * s010 - q0w * q0y * s020 -
						qperpw * qperpy * s020 - qperpx * qperpz * s020 + s100 -
						q0y * q0y * s100 - q0z * q0z * s100 - qperpy * qperpy * s100 -
						qperpz * qperpz * s100 - q0w * q0z * s110 +
						qperpx * qperpy * s110 - qperpw * qperpz * s110 +
						q0w * q0y * s120 + qperpw * qperpy * s120 +
						qperpx * qperpz * s120 +
						q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
				(-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
								s001 +
						q0w * q0z * s011 - qperpx * qperpy * s011 +
						qperpw * qperpz * s011 - q0w * q0y * s021 -
						qperpw * qperpy * s021 - qperpx * qperpz * s021 + s101 -
						q0y * q0y * s101 - q0z * q0z * s101 - qperpy * qperpy * s101 -
						qperpz * qperpz * s101 - q0w * q0z * s111 +
						qperpx * qperpy * s111 - qperpw * qperpz * s111 +
						q0w * q0y * s121 + qperpw * qperpy * s121 +
						qperpx * qperpz * s121 +
						q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
				(-1 + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) *
								s002 +
						q0w * q0z * s012 - qperpx * qperpy * s012 +
						qperpw * qperpz * s012 - q0w * q0y * s022 -
						qperpw * qperpy * s022 - qperpx * qperpz * s022 + s102 -
						q0y * q0y * s102 - q0z * q0z * s102 - qperpy * qperpy * s102 -
						qperpz * qperpz * s102 - q0w * q0z * s112 +
						qperpx * qperpy * s112 - qperpw * qperpz * s112 +
						q0w * q0y * s122 + qperpw * qperpy * s122 +
						qperpx * qperpz * s122 +
						q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122));

		c2[0] = DerivativeTerm(
				0.,
				-(qperpy * qperpy * s000) - qperpz * qperpz * s000 +
						qperpx * qperpy * s010 - qperpw * qperpz * s010 +
						qperpw * qperpy * s020 + qperpx * qperpz * s020 +
						q0y * q0y * (s000 - s100) + q0z * q0z * (s000 - s100) +
						qperpy * qperpy * s100 + qperpz * qperpz * s100 -
						qperpx * qperpy * s110 + qperpw * qperpz * s110 -
						qperpw * qperpy * s120 - qperpx * qperpz * s120 +
						2 * q0x * qperpy * s010 * theta -
						2 * q0w * qperpz * s010 * theta +
						2 * q0w * qperpy * s020 * theta +
						2 * q0x * qperpz * s020 * theta +
						q0y *
								(q0x * (-s010 + s110) + q0w * (-s020 + s120) +
										2 * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020) *
												theta) +
						q0z * (q0w * (s010 - s110) + q0x * (-s020 + s120) - 2 * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020) * theta),
				-(qperpy * qperpy * s001) - qperpz * qperpz * s001 +
						qperpx * qperpy * s011 - qperpw * qperpz * s011 +
						qperpw * qperpy * s021 + qperpx * qperpz * s021 +
						q0y * q0y * (s001 - s101) + q0z * q0z * (s001 - s101) +
						qperpy * qperpy * s101 + qperpz * qperpz * s101 -
						qperpx * qperpy * s111 + qperpw * qperpz * s111 -
						qperpw * qperpy * s121 - qperpx * qperpz * s121 +
						2 * q0x * qperpy * s011 * theta -
						2 * q0w * qperpz * s011 * theta +
						2 * q0w * qperpy * s021 * theta +
						2 * q0x * qperpz * s021 * theta +
						q0y *
								(q0x * (-s011 + s111) + q0w * (-s021 + s121) +
										2 * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021) *
												theta) +
						q0z * (q0w * (s011 - s111) + q0x * (-s021 + s121) - 2 * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021) * theta),
				-(qperpy * qperpy * s002) - qperpz * qperpz * s002 +
						qperpx * qperpy * s012 - qperpw * qperpz * s012 +
						qperpw * qperpy * s022 + qperpx * qperpz * s022 +
						q0y * q0y * (s002 - s102) + q0z * q0z * (s002 - s102) +
						qperpy * qperpy * s102 + qperpz * qperpz * s102 -
						qperpx * qperpy * s112 + qperpw * qperpz * s112 -
						qperpw * qperpy * s122 - qperpx * qperpz * s122 +
						2 * q0x * qperpy * s012 * theta -
						2 * q0w * qperpz * s012 * theta +
						2 * q0w * qperpy * s022 * theta +
						2 * q0x * qperpz * s022 * theta +
						q0y *
								(q0x * (-s012 + s112) + q0w * (-s022 + s122) +
										2 * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022) *
												theta) +
						q0z * (q0w * (s012 - s112) + q0x * (-s022 + s122) - 2 * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022) * theta));

		c3[0] = DerivativeTerm(
				0.,
				-2 * (q0x * qperpy * s010 - q0w * qperpz * s010 + q0w * qperpy * s020 + q0x * qperpz * s020 - q0x * qperpy * s110 + q0w * qperpz * s110 - q0w * qperpy * s120 - q0x * qperpz * s120 + q0y * (-2 * qperpy * s000 + qperpx * s010 + qperpw * s020 + 2 * qperpy * s100 - qperpx * s110 - qperpw * s120) + q0z * (-2 * qperpz * s000 - qperpw * s010 + qperpx * s020 + 2 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
						theta,
				-2 * (q0x * qperpy * s011 - q0w * qperpz * s011 + q0w * qperpy * s021 + q0x * qperpz * s021 - q0x * qperpy * s111 + q0w * qperpz * s111 - q0w * qperpy * s121 - q0x * qperpz * s121 + q0y * (-2 * qperpy * s001 + qperpx * s011 + qperpw * s021 + 2 * qperpy * s101 - qperpx * s111 - qperpw * s121) + q0z * (-2 * qperpz * s001 - qperpw * s011 + qperpx * s021 + 2 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
						theta,
				-2 * (q0x * qperpy * s012 - q0w * qperpz * s012 + q0w * qperpy * s022 + q0x * qperpz * s022 - q0x * qperpy * s112 + q0w * qperpz * s112 - q0w * qperpy * s122 - q0x * qperpz * s122 + q0y * (-2 * qperpy * s002 + qperpx * s012 + qperpw * s022 + 2 * qperpy * s102 - qperpx * s112 - qperpw * s122) + q0z * (-2 * qperpz * s002 - qperpw * s012 + qperpx * s022 + 2 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
						theta);

		c4[0] = DerivativeTerm(
				0.,
				-(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 -
						q0x * qperpz * s020 + q0x * qperpy * s110 -
						q0w * qperpz * s110 + q0w * qperpy * s120 +
						q0x * qperpz * s120 + 2 * q0y * q0y * s000 * theta +
						2 * q0z * q0z * s000 * theta -
						2 * qperpy * qperpy * s000 * theta -
						2 * qperpz * qperpz * s000 * theta +
						2 * qperpx * qperpy * s010 * theta -
						2 * qperpw * qperpz * s010 * theta +
						2 * qperpw * qperpy * s020 * theta +
						2 * qperpx * qperpz * s020 * theta +
						q0y * (-(qperpx * s010) - qperpw * s020 + 2 * qperpy * (s000 - s100) + qperpx * s110 + qperpw * s120 - 2 * q0x * s010 * theta - 2 * q0w * s020 * theta) +
						q0z * (2 * qperpz * s000 + qperpw * s010 - qperpx * s020 - 2 * qperpz * s100 - qperpw * s110 + qperpx * s120 + 2 * q0w * s010 * theta - 2 * q0x * s020 * theta),
				-(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 -
						q0x * qperpz * s021 + q0x * qperpy * s111 -
						q0w * qperpz * s111 + q0w * qperpy * s121 +
						q0x * qperpz * s121 + 2 * q0y * q0y * s001 * theta +
						2 * q0z * q0z * s001 * theta -
						2 * qperpy * qperpy * s001 * theta -
						2 * qperpz * qperpz * s001 * theta +
						2 * qperpx * qperpy * s011 * theta -
						2 * qperpw * qperpz * s011 * theta +
						2 * qperpw * qperpy * s021 * theta +
						2 * qperpx * qperpz * s021 * theta +
						q0y * (-(qperpx * s011) - qperpw * s021 + 2 * qperpy * (s001 - s101) + qperpx * s111 + qperpw * s121 - 2 * q0x * s011 * theta - 2 * q0w * s021 * theta) +
						q0z * (2 * qperpz * s001 + qperpw * s011 - qperpx * s021 - 2 * qperpz * s101 - qperpw * s111 + qperpx * s121 + 2 * q0w * s011 * theta - 2 * q0x * s021 * theta),
				-(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 -
						q0x * qperpz * s022 + q0x * qperpy * s112 -
						q0w * qperpz * s112 + q0w * qperpy * s122 +
						q0x * qperpz * s122 + 2 * q0y * q0y * s002 * theta +
						2 * q0z * q0z * s002 * theta -
						2 * qperpy * qperpy * s002 * theta -
						2 * qperpz * qperpz * s002 * theta +
						2 * qperpx * qperpy * s012 * theta -
						2 * qperpw * qperpz * s012 * theta +
						2 * qperpw * qperpy * s022 * theta +
						2 * qperpx * qperpz * s022 * theta +
						q0y * (-(qperpx * s012) - qperpw * s022 + 2 * qperpy * (s002 - s102) + qperpx * s112 + qperpw * s122 - 2 * q0x * s012 * theta - 2 * q0w * s022 * theta) +
						q0z * (2 * qperpz * s002 + qperpw * s012 - qperpx * s022 - 2 * qperpz * s102 - qperpw * s112 + qperpx * s122 + 2 * q0w * s012 * theta - 2 * q0x * s022 * theta));

		c5[0] = DerivativeTerm(
				0.,
				2 * (qperpy * qperpy * s000 + qperpz * qperpz * s000 - qperpx * qperpy * s010 + qperpw * qperpz * s010 - qperpw * qperpy * s020 - qperpx * qperpz * s020 - qperpy * qperpy * s100 - qperpz * qperpz * s100 + q0y * q0y * (-s000 + s100) + q0z * q0z * (-s000 + s100) + qperpx * qperpy * s110 - qperpw * qperpz * s110 + q0y * (q0x * (s010 - s110) + q0w * (s020 - s120)) + qperpw * qperpy * s120 + qperpx * qperpz * s120 + q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120)) *
						theta,
				2 * (qperpy * qperpy * s001 + qperpz * qperpz * s001 - qperpx * qperpy * s011 + qperpw * qperpz * s011 - qperpw * qperpy * s021 - qperpx * qperpz * s021 - qperpy * qperpy * s101 - qperpz * qperpz * s101 + q0y * q0y * (-s001 + s101) + q0z * q0z * (-s001 + s101) + qperpx * qperpy * s111 - qperpw * qperpz * s111 + q0y * (q0x * (s011 - s111) + q0w * (s021 - s121)) + qperpw * qperpy * s121 + qperpx * qperpz * s121 + q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121)) *
						theta,
				2 * (qperpy * qperpy * s002 + qperpz * qperpz * s002 - qperpx * qperpy * s012 + qperpw * qperpz * s012 - qperpw * qperpy * s022 - qperpx * qperpz * s022 - qperpy * qperpy * s102 - qperpz * qperpz * s102 + q0y * q0y * (-s002 + s102) + q0z * q0z * (-s002 + s102) + qperpx * qperpy * s112 - qperpw * qperpz * s112 + q0y * (q0x * (s012 - s112) + q0w * (s022 - s122)) + qperpw * qperpy * s122 + qperpx * qperpz * s122 + q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122)) *
						theta);

		c1[1] = DerivativeTerm(
				-t0y + t1y,
				-(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 +
						q0z * q0z * s010 + qperpx * qperpx * s010 +
						qperpz * qperpz * s010 - q0y * q0z * s020 +
						qperpw * qperpx * s020 - qperpy * qperpz * s020 +
						qperpx * qperpy * s100 + qperpw * qperpz * s100 +
						q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) + s110 -
						q0z * q0z * s110 - qperpx * qperpx * s110 -
						qperpz * qperpz * s110 +
						q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) +
						q0y * q0z * s120 - qperpw * qperpx * s120 +
						qperpy * qperpz * s120,
				-(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 +
						q0z * q0z * s011 + qperpx * qperpx * s011 +
						qperpz * qperpz * s011 - q0y * q0z * s021 +
						qperpw * qperpx * s021 - qperpy * qperpz * s021 +
						qperpx * qperpy * s101 + qperpw * qperpz * s101 +
						q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) + s111 -
						q0z * q0z * s111 - qperpx * qperpx * s111 -
						qperpz * qperpz * s111 +
						q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) +
						q0y * q0z * s121 - qperpw * qperpx * s121 +
						qperpy * qperpz * s121,
				-(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 +
						q0z * q0z * s012 + qperpx * qperpx * s012 +
						qperpz * qperpz * s012 - q0y * q0z * s022 +
						qperpw * qperpx * s022 - qperpy * qperpz * s022 +
						qperpx * qperpy * s102 + qperpw * qperpz * s102 +
						q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) + s112 -
						q0z * q0z * s112 - qperpx * qperpx * s112 -
						qperpz * qperpz * s112 +
						q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) +
						q0y * q0z * s122 - qperpw * qperpx * s122 +
						qperpy * qperpz * s122);

		c2[1] = DerivativeTerm(
				0.,
				qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 -
						qperpx * qperpx * s010 - qperpz * qperpz * s010 -
						q0y * q0z * s020 - qperpw * qperpx * s020 +
						qperpy * qperpz * s020 - qperpx * qperpy * s100 -
						qperpw * qperpz * s100 + q0x * q0x * (s010 - s110) -
						q0z * q0z * s110 + qperpx * qperpx * s110 +
						qperpz * qperpz * s110 + q0y * q0z * s120 +
						qperpw * qperpx * s120 - qperpy * qperpz * s120 +
						2 * q0z * qperpw * s000 * theta +
						2 * q0y * qperpx * s000 * theta -
						4 * q0z * qperpz * s010 * theta +
						2 * q0z * qperpy * s020 * theta +
						2 * q0y * qperpz * s020 * theta +
						q0x * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 + 2 * qperpy * s000 * theta - 4 * qperpx * s010 * theta - 2 * qperpw * s020 * theta) +
						q0w * (-(q0z * s000) + q0z * s100 + 2 * qperpz * s000 * theta - 2 * qperpx * s020 * theta),
				qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 -
						qperpx * qperpx * s011 - qperpz * qperpz * s011 -
						q0y * q0z * s021 - qperpw * qperpx * s021 +
						qperpy * qperpz * s021 - qperpx * qperpy * s101 -
						qperpw * qperpz * s101 + q0x * q0x * (s011 - s111) -
						q0z * q0z * s111 + qperpx * qperpx * s111 +
						qperpz * qperpz * s111 + q0y * q0z * s121 +
						qperpw * qperpx * s121 - qperpy * qperpz * s121 +
						2 * q0z * qperpw * s001 * theta +
						2 * q0y * qperpx * s001 * theta -
						4 * q0z * qperpz * s011 * theta +
						2 * q0z * qperpy * s021 * theta +
						2 * q0y * qperpz * s021 * theta +
						q0x * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 + 2 * qperpy * s001 * theta - 4 * qperpx * s011 * theta - 2 * qperpw * s021 * theta) +
						q0w * (-(q0z * s001) + q0z * s101 + 2 * qperpz * s001 * theta - 2 * qperpx * s021 * theta),
				qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 -
						qperpx * qperpx * s012 - qperpz * qperpz * s012 -
						q0y * q0z * s022 - qperpw * qperpx * s022 +
						qperpy * qperpz * s022 - qperpx * qperpy * s102 -
						qperpw * qperpz * s102 + q0x * q0x * (s012 - s112) -
						q0z * q0z * s112 + qperpx * qperpx * s112 +
						qperpz * qperpz * s112 + q0y * q0z * s122 +
						qperpw * qperpx * s122 - qperpy * qperpz * s122 +
						2 * q0z * qperpw * s002 * theta +
						2 * q0y * qperpx * s002 * theta -
						4 * q0z * qperpz * s012 * theta +
						2 * q0z * qperpy * s022 * theta +
						2 * q0y * qperpz * s022 * theta +
						q0x * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 + 2 * qperpy * s002 * theta - 4 * qperpx * s012 * theta - 2 * qperpw * s022 * theta) +
						q0w * (-(q0z * s002) + q0z * s102 + 2 * qperpz * s002 * theta - 2 * qperpx * s022 * theta));

		c3[1] = DerivativeTerm(
				0., 2 * (-(q0x * qperpy * s000) - q0w * qperpz * s000 + 2 * q0x * qperpx * s010 + q0x * qperpw * s020 + q0w * qperpx * s020 + q0x * qperpy * s100 + q0w * qperpz * s100 - 2 * q0x * qperpx * s110 - q0x * qperpw * s120 - q0w * qperpx * s120 + q0z * (2 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) - 2 * qperpz * s110 + qperpy * s120) + q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120)) * theta,
				2 * (-(q0x * qperpy * s001) - q0w * qperpz * s001 + 2 * q0x * qperpx * s011 + q0x * qperpw * s021 + q0w * qperpx * s021 + q0x * qperpy * s101 + q0w * qperpz * s101 - 2 * q0x * qperpx * s111 - q0x * qperpw * s121 - q0w * qperpx * s121 + q0z * (2 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) - 2 * qperpz * s111 + qperpy * s121) + q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121)) *
						theta,
				2 * (-(q0x * qperpy * s002) - q0w * qperpz * s002 + 2 * q0x * qperpx * s012 + q0x * qperpw * s022 + q0w * qperpx * s022 + q0x * qperpy * s102 + q0w * qperpz * s102 - 2 * q0x * qperpx * s112 - q0x * qperpw * s122 - q0w * qperpx * s122 + q0z * (2 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) - 2 * qperpz * s112 + qperpy * s122) + q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122)) *
						theta);

		c4[1] = DerivativeTerm(
				0.,
				-(q0x * qperpy * s000) - q0w * qperpz * s000 +
						2 * q0x * qperpx * s010 + q0x * qperpw * s020 +
						q0w * qperpx * s020 + q0x * qperpy * s100 +
						q0w * qperpz * s100 - 2 * q0x * qperpx * s110 -
						q0x * qperpw * s120 - q0w * qperpx * s120 +
						2 * qperpx * qperpy * s000 * theta +
						2 * qperpw * qperpz * s000 * theta +
						2 * q0x * q0x * s010 * theta + 2 * q0z * q0z * s010 * theta -
						2 * qperpx * qperpx * s010 * theta -
						2 * qperpz * qperpz * s010 * theta +
						2 * q0w * q0x * s020 * theta -
						2 * qperpw * qperpx * s020 * theta +
						2 * qperpy * qperpz * s020 * theta +
						q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120 - 2 * q0x * s000 * theta) +
						q0z * (2 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) - 2 * qperpz * s110 + qperpy * s120 - 2 * q0w * s000 * theta - 2 * q0y * s020 * theta),
				-(q0x * qperpy * s001) - q0w * qperpz * s001 +
						2 * q0x * qperpx * s011 + q0x * qperpw * s021 +
						q0w * qperpx * s021 + q0x * qperpy * s101 +
						q0w * qperpz * s101 - 2 * q0x * qperpx * s111 -
						q0x * qperpw * s121 - q0w * qperpx * s121 +
						2 * qperpx * qperpy * s001 * theta +
						2 * qperpw * qperpz * s001 * theta +
						2 * q0x * q0x * s011 * theta + 2 * q0z * q0z * s011 * theta -
						2 * qperpx * qperpx * s011 * theta -
						2 * qperpz * qperpz * s011 * theta +
						2 * q0w * q0x * s021 * theta -
						2 * qperpw * qperpx * s021 * theta +
						2 * qperpy * qperpz * s021 * theta +
						q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121 - 2 * q0x * s001 * theta) +
						q0z * (2 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) - 2 * qperpz * s111 + qperpy * s121 - 2 * q0w * s001 * theta - 2 * q0y * s021 * theta),
				-(q0x * qperpy * s002) - q0w * qperpz * s002 +
						2 * q0x * qperpx * s012 + q0x * qperpw * s022 +
						q0w * qperpx * s022 + q0x * qperpy * s102 +
						q0w * qperpz * s102 - 2 * q0x * qperpx * s112 -
						q0x * qperpw * s122 - q0w * qperpx * s122 +
						2 * qperpx * qperpy * s002 * theta +
						2 * qperpw * qperpz * s002 * theta +
						2 * q0x * q0x * s012 * theta + 2 * q0z * q0z * s012 * theta -
						2 * qperpx * qperpx * s012 * theta -
						2 * qperpz * qperpz * s012 * theta +
						2 * q0w * q0x * s022 * theta -
						2 * qperpw * qperpx * s022 * theta +
						2 * qperpy * qperpz * s022 * theta +
						q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122 - 2 * q0x * s002 * theta) +
						q0z * (2 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) - 2 * qperpz * s112 + qperpy * s122 - 2 * q0w * s002 * theta - 2 * q0y * s022 * theta));

		c5[1] = DerivativeTerm(
				0., -2 * (qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010 - qperpx * qperpx * s010 - qperpz * qperpz * s010 - q0y * q0z * s020 - qperpw * qperpx * s020 + qperpy * qperpz * s020 - qperpx * qperpy * s100 - qperpw * qperpz * s100 + q0w * q0z * (-s000 + s100) + q0x * q0x * (s010 - s110) - q0z * q0z * s110 + qperpx * qperpx * s110 + qperpz * qperpz * s110 + q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120)) + q0y * q0z * s120 + qperpw * qperpx * s120 - qperpy * qperpz * s120) * theta,
				-2 * (qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011 - qperpx * qperpx * s011 - qperpz * qperpz * s011 - q0y * q0z * s021 - qperpw * qperpx * s021 + qperpy * qperpz * s021 - qperpx * qperpy * s101 - qperpw * qperpz * s101 + q0w * q0z * (-s001 + s101) + q0x * q0x * (s011 - s111) - q0z * q0z * s111 + qperpx * qperpx * s111 + qperpz * qperpz * s111 + q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121)) + q0y * q0z * s121 + qperpw * qperpx * s121 - qperpy * qperpz * s121) *
						theta,
				-2 * (qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012 - qperpx * qperpx * s012 - qperpz * qperpz * s012 - q0y * q0z * s022 - qperpw * qperpx * s022 + qperpy * qperpz * s022 - qperpx * qperpy * s102 - qperpw * qperpz * s102 + q0w * q0z * (-s002 + s102) + q0x * q0x * (s012 - s112) - q0z * q0z * s112 + qperpx * qperpx * s112 + qperpz * qperpz * s112 + q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122)) + q0y * q0z * s122 + qperpw * qperpx * s122 - qperpy * qperpz * s122) *
						theta);

		c1[2] = DerivativeTerm(
				-t0z + t1z, (qperpw * qperpy * s000 - qperpx * qperpz * s000 - q0y * q0z * s010 - qperpw * qperpx * s010 - qperpy * qperpz * s010 - s020 + q0y * q0y * s020 + qperpx * qperpx * s020 + qperpy * qperpy * s020 - qperpw * qperpy * s100 + qperpx * qperpz * s100 + q0x * q0z * (-s000 + s100) + q0y * q0z * s110 + qperpw * qperpx * s110 + qperpy * qperpz * s110 + q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110)) + q0x * q0x * (s020 - s120) + s120 - q0y * q0y * s120 - qperpx * qperpx * s120 - qperpy * qperpy * s120),
				(qperpw * qperpy * s001 - qperpx * qperpz * s001 -
						q0y * q0z * s011 - qperpw * qperpx * s011 -
						qperpy * qperpz * s011 - s021 + q0y * q0y * s021 +
						qperpx * qperpx * s021 + qperpy * qperpy * s021 -
						qperpw * qperpy * s101 + qperpx * qperpz * s101 +
						q0x * q0z * (-s001 + s101) + q0y * q0z * s111 +
						qperpw * qperpx * s111 + qperpy * qperpz * s111 +
						q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111)) +
						q0x * q0x * (s021 - s121) + s121 - q0y * q0y * s121 -
						qperpx * qperpx * s121 - qperpy * qperpy * s121),
				(qperpw * qperpy * s002 - qperpx * qperpz * s002 -
						q0y * q0z * s012 - qperpw * qperpx * s012 -
						qperpy * qperpz * s012 - s022 + q0y * q0y * s022 +
						qperpx * qperpx * s022 + qperpy * qperpy * s022 -
						qperpw * qperpy * s102 + qperpx * qperpz * s102 +
						q0x * q0z * (-s002 + s102) + q0y * q0z * s112 +
						qperpw * qperpx * s112 + qperpy * qperpz * s112 +
						q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112)) +
						q0x * q0x * (s022 - s122) + s122 - q0y * q0y * s122 -
						qperpx * qperpx * s122 - qperpy * qperpy * s122));

		c2[2] = DerivativeTerm(
				0.,
				(q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 +
						qperpx * qperpz * s000 - q0w * q0x * s010 - q0y * q0z * s010 +
						qperpw * qperpx * s010 + qperpy * qperpz * s010 +
						q0x * q0x * s020 + q0y * q0y * s020 - qperpx * qperpx * s020 -
						qperpy * qperpy * s020 - q0w * q0y * s100 + q0x * q0z * s100 +
						qperpw * qperpy * s100 - qperpx * qperpz * s100 +
						q0w * q0x * s110 + q0y * q0z * s110 - qperpw * qperpx * s110 -
						qperpy * qperpz * s110 - q0x * q0x * s120 - q0y * q0y * s120 +
						qperpx * qperpx * s120 + qperpy * qperpy * s120 -
						2 * q0y * qperpw * s000 * theta + 2 * q0z * qperpx * s000 * theta -
						2 * q0w * qperpy * s000 * theta + 2 * q0x * qperpz * s000 * theta +
						2 * q0x * qperpw * s010 * theta + 2 * q0w * qperpx * s010 * theta +
						2 * q0z * qperpy * s010 * theta + 2 * q0y * qperpz * s010 * theta -
						4 * q0x * qperpx * s020 * theta - 4 * q0y * qperpy * s020 * theta),
				(q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 +
						qperpx * qperpz * s001 - q0w * q0x * s011 - q0y * q0z * s011 +
						qperpw * qperpx * s011 + qperpy * qperpz * s011 +
						q0x * q0x * s021 + q0y * q0y * s021 - qperpx * qperpx * s021 -
						qperpy * qperpy * s021 - q0w * q0y * s101 + q0x * q0z * s101 +
						qperpw * qperpy * s101 - qperpx * qperpz * s101 +
						q0w * q0x * s111 + q0y * q0z * s111 - qperpw * qperpx * s111 -
						qperpy * qperpz * s111 - q0x * q0x * s121 - q0y * q0y * s121 +
						qperpx * qperpx * s121 + qperpy * qperpy * s121 -
						2 * q0y * qperpw * s001 * theta + 2 * q0z * qperpx * s001 * theta -
						2 * q0w * qperpy * s001 * theta + 2 * q0x * qperpz * s001 * theta +
						2 * q0x * qperpw * s011 * theta + 2 * q0w * qperpx * s011 * theta +
						2 * q0z * qperpy * s011 * theta + 2 * q0y * qperpz * s011 * theta -
						4 * q0x * qperpx * s021 * theta - 4 * q0y * qperpy * s021 * theta),
				(q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 +
						qperpx * qperpz * s002 - q0w * q0x * s012 - q0y * q0z * s012 +
						qperpw * qperpx * s012 + qperpy * qperpz * s012 +
						q0x * q0x * s022 + q0y * q0y * s022 - qperpx * qperpx * s022 -
						qperpy * qperpy * s022 - q0w * q0y * s102 + q0x * q0z * s102 +
						qperpw * qperpy * s102 - qperpx * qperpz * s102 +
						q0w * q0x * s112 + q0y * q0z * s112 - qperpw * qperpx * s112 -
						qperpy * qperpz * s112 - q0x * q0x * s122 - q0y * q0y * s122 +
						qperpx * qperpx * s122 + qperpy * qperpy * s122 -
						2 * q0y * qperpw * s002 * theta + 2 * q0z * qperpx * s002 * theta -
						2 * q0w * qperpy * s002 * theta + 2 * q0x * qperpz * s002 * theta +
						2 * q0x * qperpw * s012 * theta + 2 * q0w * qperpx * s012 * theta +
						2 * q0z * qperpy * s012 * theta + 2 * q0y * qperpz * s012 * theta -
						4 * q0x * qperpx * s022 * theta -
						4 * q0y * qperpy * s022 * theta));

		c3[2] = DerivativeTerm(
				0., -2 * (-(q0w * qperpy * s000) + q0x * qperpz * s000 + q0x * qperpw * s010 + q0w * qperpx * s010 - 2 * q0x * qperpx * s020 + q0w * qperpy * s100 - q0x * qperpz * s100 - q0x * qperpw * s110 - q0w * qperpx * s110 + q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 - qperpy * s110) + 2 * q0x * qperpx * s120 + q0y * (qperpz * s010 - 2 * qperpy * s020 + qperpw * (-s000 + s100) - qperpz * s110 + 2 * qperpy * s120)) * theta,
				-2 * (-(q0w * qperpy * s001) + q0x * qperpz * s001 + q0x * qperpw * s011 + q0w * qperpx * s011 - 2 * q0x * qperpx * s021 + q0w * qperpy * s101 - q0x * qperpz * s101 - q0x * qperpw * s111 - q0w * qperpx * s111 + q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 - qperpy * s111) + 2 * q0x * qperpx * s121 + q0y * (qperpz * s011 - 2 * qperpy * s021 + qperpw * (-s001 + s101) - qperpz * s111 + 2 * qperpy * s121)) *
						theta,
				-2 * (-(q0w * qperpy * s002) + q0x * qperpz * s002 + q0x * qperpw * s012 + q0w * qperpx * s012 - 2 * q0x * qperpx * s022 + q0w * qperpy * s102 - q0x * qperpz * s102 - q0x * qperpw * s112 - q0w * qperpx * s112 + q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 - qperpy * s112) + 2 * q0x * qperpx * s122 + q0y * (qperpz * s012 - 2 * qperpy * s022 + qperpw * (-s002 + s102) - qperpz * s112 + 2 * qperpy * s122)) *
						theta);

		c4[2] = DerivativeTerm(
				0.,
				q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 -
						q0w * qperpx * s010 + 2 * q0x * qperpx * s020 -
						q0w * qperpy * s100 + q0x * qperpz * s100 +
						q0x * qperpw * s110 + q0w * qperpx * s110 -
						2 * q0x * qperpx * s120 - 2 * qperpw * qperpy * s000 * theta +
						2 * qperpx * qperpz * s000 * theta -
						2 * q0w * q0x * s010 * theta +
						2 * qperpw * qperpx * s010 * theta +
						2 * qperpy * qperpz * s010 * theta +
						2 * q0x * q0x * s020 * theta + 2 * q0y * q0y * s020 * theta -
						2 * qperpx * qperpx * s020 * theta -
						2 * qperpy * qperpy * s020 * theta +
						q0z * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 + qperpy * s110 - 2 * q0x * s000 * theta) +
						q0y * (-(qperpz * s010) + 2 * qperpy * s020 + qperpw * (s000 - s100) + qperpz * s110 - 2 * qperpy * s120 + 2 * q0w * s000 * theta - 2 * q0z * s010 * theta),
				q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 -
						q0w * qperpx * s011 + 2 * q0x * qperpx * s021 -
						q0w * qperpy * s101 + q0x * qperpz * s101 +
						q0x * qperpw * s111 + q0w * qperpx * s111 -
						2 * q0x * qperpx * s121 - 2 * qperpw * qperpy * s001 * theta +
						2 * qperpx * qperpz * s001 * theta -
						2 * q0w * q0x * s011 * theta +
						2 * qperpw * qperpx * s011 * theta +
						2 * qperpy * qperpz * s011 * theta +
						2 * q0x * q0x * s021 * theta + 2 * q0y * q0y * s021 * theta -
						2 * qperpx * qperpx * s021 * theta -
						2 * qperpy * qperpy * s021 * theta +
						q0z * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 + qperpy * s111 - 2 * q0x * s001 * theta) +
						q0y * (-(qperpz * s011) + 2 * qperpy * s021 + qperpw * (s001 - s101) + qperpz * s111 - 2 * qperpy * s121 + 2 * q0w * s001 * theta - 2 * q0z * s011 * theta),
				q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 -
						q0w * qperpx * s012 + 2 * q0x * qperpx * s022 -
						q0w * qperpy * s102 + q0x * qperpz * s102 +
						q0x * qperpw * s112 + q0w * qperpx * s112 -
						2 * q0x * qperpx * s122 - 2 * qperpw * qperpy * s002 * theta +
						2 * qperpx * qperpz * s002 * theta -
						2 * q0w * q0x * s012 * theta +
						2 * qperpw * qperpx * s012 * theta +
						2 * qperpy * qperpz * s012 * theta +
						2 * q0x * q0x * s022 * theta + 2 * q0y * q0y * s022 * theta -
						2 * qperpx * qperpx * s022 * theta -
						2 * qperpy * qperpy * s022 * theta +
						q0z * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 + qperpy * s112 - 2 * q0x * s002 * theta) +
						q0y * (-(qperpz * s012) + 2 * qperpy * s022 + qperpw * (s002 - s102) + qperpz * s112 - 2 * qperpy * s122 + 2 * q0w * s002 * theta - 2 * q0z * s012 * theta));

		c5[2] = DerivativeTerm(
				0., 2 * (qperpw * qperpy * s000 - qperpx * qperpz * s000 + q0y * q0z * s010 - qperpw * qperpx * s010 - qperpy * qperpz * s010 - q0y * q0y * s020 + qperpx * qperpx * s020 + qperpy * qperpy * s020 + q0x * q0z * (s000 - s100) - qperpw * qperpy * s100 + qperpx * qperpz * s100 + q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110)) - q0y * q0z * s110 + qperpw * qperpx * s110 + qperpy * qperpz * s110 + q0y * q0y * s120 - qperpx * qperpx * s120 - qperpy * qperpy * s120 + q0x * q0x * (-s020 + s120)) * theta,
				2 * (qperpw * qperpy * s001 - qperpx * qperpz * s001 + q0y * q0z * s011 - qperpw * qperpx * s011 - qperpy * qperpz * s011 - q0y * q0y * s021 + qperpx * qperpx * s021 + qperpy * qperpy * s021 + q0x * q0z * (s001 - s101) - qperpw * qperpy * s101 + qperpx * qperpz * s101 + q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111)) - q0y * q0z * s111 + qperpw * qperpx * s111 + qperpy * qperpz * s111 + q0y * q0y * s121 - qperpx * qperpx * s121 - qperpy * qperpy * s121 + q0x * q0x * (-s021 + s121)) *
						theta,
				2 * (qperpw * qperpy * s002 - qperpx * qperpz * s002 + q0y * q0z * s012 - qperpw * qperpx * s012 - qperpy * qperpz * s012 - q0y * q0y * s022 + qperpx * qperpx * s022 + qperpy * qperpy * s022 + q0x * q0z * (s002 - s102) - qperpw * qperpy * s102 + qperpx * qperpz * s102 + q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112)) - q0y * q0z * s112 + qperpw * qperpx * s112 + qperpy * qperpz * s112 + q0y * q0y * s122 - qperpx * qperpx * s122 - qperpy * qperpy * s122 + q0x * q0x * (-s022 + s122)) *
						theta);
	}
}

void AnimatedTransform::decompose(const Matrix4& m, Vector3f* T, Quaternion* Rquat, Matrix4* S) {
	// extract translation
	T->x = m.m[0][3];
	T->y = m.m[1][3];
	T->z = m.m[2][3];

	// assume affine, so only need upper 3x3 matrix
	Matrix4 M = m;
	for (int i = 0; i < 3; i++) {
		M.m[i][3] = M.m[3][i] = 0.0f;
	}
	M.m[3][3] = 1.0f;

	// extract rotation by averaging M with its inverse transpose until converge
	float norm;
	int count = 0;
	Matrix4 R = M;
	do {
		Matrix4 Rnext;
		Matrix4 Rit = inverse(transpose(R));
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Rnext.m[i][j] = 0.5 * (R.m[i][j] + Rit.m[i][j]);
			}
		}

		norm = 0;
		for (int i = 0; i < 3; i++) {
			float n = std::abs(R.m[i][0] - Rnext.m[i][0]) +
					std::abs(R.m[i][1] - Rnext.m[i][1]) +
					std::abs(R.m[i][2] - Rnext.m[i][2]);
			norm = std::max(norm, n);
		}

		R = Rnext;
	} while (++count < 100 && norm > 0.0001f);
	*Rquat = Quaternion(R);

	// extract scale using rotation and original M
	*S = Matrix4::matmul(inverse(R), M);
}

void AnimatedTransform::interpolate(float time, Transform* t) const {
	if (!is_animated || time <= start_time) {
		*t = *start_transform;
		return;
	}
	if (time >= end_time) {
		*t = *end_transform;
		return;
	}

	float dt = (time - start_time) / (end_time - start_time);
	Vector3f trans = (1 - dt) * T[0] + dt * T[1];
	Quaternion rotate = slerp(dt, R[0], R[1]);
	Matrix4 scale;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			scale.m[i][j] = lerp(dt, S[0].m[i][j], S[1].m[i][j]);
		}
	}

	*t = translate(trans) * rotate.toTransform() * Transform(scale);
}

Ray AnimatedTransform::operator()(const Ray& r) const {
	if (!is_animated || r.time <= start_time) {
		return (*start_transform)(r);
	} else if (r.time >= end_time) {
		return (*end_transform)(r);
	} else {
		Transform t;
		interpolate(r.time, &t);
		return t(r);
	}
}

RayDifferential AnimatedTransform::operator()(const RayDifferential& r) const {
	if (!is_animated || r.time <= start_time) {
		return (*start_transform)(r);
	} else if (r.time >= end_time) {
		return (*end_transform)(r);
	} else {
		Transform t;
		interpolate(r.time, &t);
		return t(r);
	}
}

Point3f AnimatedTransform::operator()(float time, const Point3f& p) const {
	if (!is_animated || time <= start_time) {
		return (*start_transform)(p);
	} else if (time >= end_time) {
		return (*end_transform)(p);
	} else {
		Transform t;
		interpolate(time, &t);
		return t(p);
	}
}

Vector3f AnimatedTransform::operator()(float time, const Vector3f& v) const {
	if (!is_animated || time <= start_time) {
		return (*start_transform)(v);
	} else if (time >= end_time) {
		return (*end_transform)(v);
	} else {
		Transform t;
		interpolate(time, &t);
		return t(v);
	}
}

Bounds3f AnimatedTransform::motionBounds(const Bounds3f& b) const {
	if (!is_animated) {
		return (*start_transform)(b);
	}
	if (!has_rotation) {
		return bUnion((*start_transform)(b), (*end_transform)(b));
	}

	Bounds3f bounds;
	for (int corner = 0; corner < 8; corner++) {
		bounds = bUnion(bounds, boundPointMotion(b.corner(corner)));
	}
	return bounds;
}

Bounds3f AnimatedTransform::boundPointMotion(const Point3f& p) const {
	Bounds3f bounds((*start_transform)(p), (*end_transform)(p));
	float cos_theta = dot(R[0], R[1]);
	float theta = acos(clamp(cos_theta, -1, 1));
	for (int c = 0; c < 3; c++) {
		float zeros[4];
		int n_zeros = 0;
		intervalFindZeros(c1[c].eval(p), c2[c].eval(p), c3[c].eval(p), c4[c].eval(p), c5[c].eval(p),
				theta, Interval(0, 1), zeros, &n_zeros);

		// expand bounds for any derivative zeros found
		for (int i = 0; i < n_zeros; i++) {
			Point3f pz = (*this)(lerp(zeros[i], start_time, end_time), p);
			bounds = bUnion(bounds, pz);
		}
	}
	return bounds;
}

} // namespace civet