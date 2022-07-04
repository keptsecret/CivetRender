#include <shapes/sphere.h>

#include <core/interaction.h>	///< TODO: this might have to be removed

namespace civet {

CIVET_CPU_GPU
Bounds3f Sphere::objectBound() const {
	return Bounds3f(Point3f(-radius, -radius, z_min), Point3f(radius, radius, z_max));
}

CIVET_CPU_GPU
bool Sphere::intersect(const Ray& ray, float* t_hit, SurfaceInteraction* isect, bool test_alpha_texture) const {
	///< Maybe will implement the floating point bound error version, will see results first
	// Vector3f o_err, d_err;
	Ray r = (*world_to_object)(ray);

	// Compute quadratic sphere coefficients
	// TODO: implement and use EFloat class here, maybe
	float ox = r.o.x, oy = r.o.y, oz = r.o.z;
	float dx = r.d.x, dy = r.d.y, dz = r.d.z;

	float a = dx * dx + dy * dy + dz * dz;
	float b = 2 * (dx * ox + dy * oy + dz * oz);
	float c = ox * ox + oy * oy + oz * oz - radius * radius;

	// Solve quadratic equation for t values
	float t0, t1;
	if (!quadratic(a, b, c, t0, t1)) {
		return false;
	}

	// Check nearest intersection
	if (t0 > r.t_max || t1 <= 0) {
		return false;
	}
	float t_shape_hit = t0;
	if (t_shape_hit <= 0) {
		t_shape_hit = t1;
		if (t_shape_hit > r.t_max) {
			return false;
		}
	}

	// compute hit position and phi
	float phi;
	Point3f p_hit;
	p_hit = ray(t_shape_hit);
	p_hit *= radius / distance(p_hit, Point3f(0, 0, 0)); ///< account for precision limitations
	if (p_hit.x == 0 && p_hit.y == 0) {
		p_hit.x = 1e-5f * radius;
	}
	phi = std::atan2(p_hit.y, p_hit.x);
	if (phi < 0) {
		phi += 2 * Pi;
	}

	// test intersection against clipping params, z_min/max and phi
	if ((z_min > -radius && p_hit.z < z_min) || (z_max < radius && p_hit.z > z_max) || phi > phi_max) {
		if (t_shape_hit == t1) {
			return false;
		}
		if (t1 > r.t_max) {
			return false;
		}
		t_shape_hit = t1;

		// compute again for t1
		p_hit = ray(t_shape_hit);
		p_hit *= radius / distance(p_hit, Point3f(0, 0, 0)); ///< account for precision limitations
		if (p_hit.x == 0 && p_hit.y == 0) {
			p_hit.x = 1e-5f * radius;
		}
		phi = std::atan2(p_hit.y, p_hit.x);
		if (phi < 0) {
			phi += 2 * Pi;
		}
		if ((z_min > -radius && p_hit.z < z_min) || (z_max < radius && p_hit.z > z_max) || phi > phi_max) {
			return false;
		}
	}

	// confirmed hit, calculate parametric representation of sphere
	float u = phi / phi_max;
	float theta = std::acos(clamp(p_hit.z / radius, -1, 1));
	float v = (theta - theta_min) / (theta_max - theta_min);

	// compute dpdu, dpdv
	float z_radius = std::sqrt(p_hit.x * p_hit.x + p_hit.y * p_hit.y);
	float inv_z_radius = 1 / z_radius;
	float cos_phi = p_hit.x * inv_z_radius;
	float sin_phi = p_hit.y * inv_z_radius;
	Vector3f dpdu(-phi_max * p_hit.y, phi_max * p_hit.x, 0);
	Vector3f dpdv = (theta_max - theta_min) * Vector3f(p_hit.z * cos_phi, p_hit.z * sin_phi, -radius * std::sin(theta));

	// compute dndu, dndv
	Vector3f d2pduu = -phi_max * phi_max * Vector3f(p_hit.x, p_hit.y, 0);
	Vector3f d2pduv = (theta_max - theta_min) * p_hit.z * phi_max * Vector3f(-sin_phi, cos_phi, 0);
	Vector3f d2pdvv = -(theta_max - theta_min) * (theta_max - theta_min) * Vector3f(p_hit);

	float E = dot(dpdu, dpdu);
	float F = dot(dpdu, dpdv);
	float G = dot(dpdv, dpdv);
	Vector3f N = normalize(cross(dpdu, dpdv));
	float e = dot(N, d2pduu);
	float f = dot(N, d2pduv);
	float g = dot(N, d2pdvv);

	float inv_EGF2 = 1 / (E * G - F * F);
	Normal3f dndu((f * F - e * G) * inv_EGF2 * dpdu + (e * F - f * E) * inv_EGF2 * dpdv);
	Normal3f dndv((g * F - f * G) * inv_EGF2 * dpdu + (f * F - g * E) * inv_EGF2 * dpdv);

	// compute error bounds for sphere intersection
	Vector3f p_error = gamma(5) * abs((Vector3f)p_hit);

	// initialize SurfaceInteraction from info
	*isect = (*object_to_world)(SurfaceInteraction(p_hit, p_error, Point2f(u, v),-r.d, dpdu, dpdv, dndu, dndv, r.time, this));
	*t_hit = t_shape_hit;
	return true;
}
bool Sphere::intersectP(const Ray& ray, bool test_alpha_texture) const {
	///< Same as intersect() but without creating interactions
	Ray r = (*world_to_object)(ray);

	// Compute quadratic sphere coefficients
	// TODO: implement and use EFloat class here, maybe
	float ox = r.o.x, oy = r.o.y, oz = r.o.z;
	float dx = r.d.x, dy = r.d.y, dz = r.d.z;

	float a = dx * dx + dy * dy + dz * dz;
	float b = 2 * (dx * ox + dy * oy + dz * oz);
	float c = ox * ox + oy * oy + oz * oz - radius * radius;

	// Solve quadratic equation for t values
	float t0, t1;
	if (!quadratic(a, b, c, t0, t1)) {
		return false;
	}

	// Check nearest intersection
	if (t0 > r.t_max || t1 <= 0) {
		return false;
	}
	float t_shape_hit = t0;
	if (t_shape_hit <= 0) {
		t_shape_hit = t1;
		if (t_shape_hit > r.t_max) {
			return false;
		}
	}

	// compute hit position and phi
	float phi;
	Point3f p_hit;
	p_hit = ray(t_shape_hit);
	p_hit *= radius / distance(p_hit, Point3f(0, 0, 0)); ///< account for precision limitations
	if (p_hit.x == 0 && p_hit.y == 0) {
		p_hit.x = 1e-5f * radius;
	}
	phi = std::atan2(p_hit.y, p_hit.x);
	if (phi < 0) {
		phi += 2 * Pi;
	}

	// test intersection against clipping params, z_min/max and phi
	if ((z_min > -radius && p_hit.z < z_min) || (z_max < radius && p_hit.z > z_max) || phi > phi_max) {
		if (t_shape_hit == t1) {
			return false;
		}
		if (t1 > r.t_max) {
			return false;
		}
		t_shape_hit = t1;

		// compute again for t1
		p_hit = ray(t_shape_hit);
		p_hit *= radius / distance(p_hit, Point3f(0, 0, 0)); ///< account for precision limitations
		if (p_hit.x == 0 && p_hit.y == 0) {
			p_hit.x = 1e-5f * radius;
		}
		phi = std::atan2(p_hit.y, p_hit.x);
		if (phi < 0) {
			phi += 2 * Pi;
		}
		if ((z_min > -radius && p_hit.z < z_min) || (z_max < radius && p_hit.z > z_max) || phi > phi_max) {
			return false;
		}
	}

	return true;
}
float Sphere::area() const {
	return phi_max * radius * (z_max - z_min);
}

} // namespace civet