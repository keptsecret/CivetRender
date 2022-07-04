#include "transform.h"

#include "interaction.h"

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
				swap(minv[irow][k], minv[icol][k]);
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
				swap(minv[k][indxr[j]], minv[k][indxc[j]]);
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
Transform lookAt(const Point3f& position, const Point3f& target, const Vector3f& up) {
	Matrix4 camera_to_world;
	camera_to_world.m[0][3] = position.x;
	camera_to_world.m[1][3] = position.y;
	camera_to_world.m[2][3] = position.z;
	camera_to_world.m[3][3] = 1;

	Vector3f dir = normalize(target - position);
	Vector3f right = normalize(cross(normalize(up), dir));
	Vector3f new_up = normalize(cross(dir, right));

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

} // namespace civet