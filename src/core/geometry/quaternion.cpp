#include <core/geometry/quaternion.h>

#include <core/geometry/transform.h>

namespace civet {

Quaternion::Quaternion(const Transform& t) {
	const Matrix4 &m = t.m;
	float trace = m.m[0][0] + m.m[1][1] + m.m[2][2];
	if (trace > 0.f) {
		float s = sqrt(trace + 1.0f);
		w = s / 2.0f;
		s = 0.5f / s;
		v.x = (m.m[2][1] - m.m[1][2]) * s;
		v.y = (m.m[0][2] - m.m[2][0]) * s;
		v.z = (m.m[1][0] - m.m[0][1]) * s;
	} else {
		const int nxt[3] = {1, 2, 0};
		float q[3];
		int i = 0;
		if (m.m[1][1] > m.m[0][0]) {
			i = 1;
		}
		if (m.m[2][2] > m.m[i][i]) {
			i = 2;
		}
		int j = nxt[i];
		int k = nxt[j];
		float s = sqrt((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + 1.0f);
		q[i] = s * 0.5f;
		if (s != 0.f) {
			s = 0.5f / s;
		}
		w = (m.m[k][j] - m.m[j][k]) * s;
		q[j] = (m.m[j][i] + m.m[i][j]) * s;
		q[k] = (m.m[k][i] + m.m[i][k]) * s;
		v.x = q[0];
		v.y = q[1];
		v.z = q[2];
	}
}

Transform Quaternion::toTransform() const {
	float xx = v.x * v.x, yy = v.y * v.y, zz = v.z * v.z;
	float xy = v.x * v.y, xz = v.x * v.z, yz = v.y * v.z;
	float wx = v.x * w, wy = v.y * w, wz = v.z * w;

	Matrix4 m;
	m.m[0][0] = 1 - 2 * (yy + zz);
	m.m[0][1] = 2 * (xy + wz);
	m.m[0][2] = 2 * (xz - wy);
	m.m[1][0] = 2 * (xy - wz);
	m.m[1][1] = 1 - 2 * (xx + zz);
	m.m[1][2] = 2 * (yz + wx);
	m.m[2][0] = 2 * (xz + wy);
	m.m[2][1] = 2 * (yz - wx);
	m.m[2][2] = 1 - 2 * (xx + yy);

	return Transform(transpose(m), m);
}

Quaternion slerp(float t, const Quaternion& q1, const Quaternion& q2) {
	float cos_theta = dot(q1, q2);
	if (cos_theta > 0.9995f) {
		return normalize((1-t)*q1 + t * q2);
	} else {
		float theta = acos(clamp(cos_theta, -1, 1));
		float thetap = theta * t;
		Quaternion qperp = normalize(q2 - q1 * cos_theta);
		return q1 * cos(thetap) + qperp * sin(thetap);
	}
}

} // namespace civet