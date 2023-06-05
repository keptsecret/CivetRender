#include <shapes/triangle.h>

#include <core/interaction.h>

namespace civet {

CIVET_CPU_GPU
Bounds3f Triangle::objectBound() const {
	const Point3f& p0 = mesh->p[v[0]];
	const Point3f& p1 = mesh->p[v[1]];
	const Point3f& p2 = mesh->p[v[2]];
	return bUnion(Bounds3f((*world_to_object)(p0), (*world_to_object)(p1)), (*world_to_object)(p2));
}

CIVET_CPU_GPU
Bounds3f Triangle::worldBound() const {
	const Point3f& p0 = mesh->p[v[0]];
	const Point3f& p1 = mesh->p[v[1]];
	const Point3f& p2 = mesh->p[v[2]];
	return bUnion(Bounds3f(p0, p1), p2);
}

CIVET_CPU_GPU
bool Triangle::intersect(const Ray& ray, float* t_hit, SurfaceInteraction* isect, bool test_alpha_texture) const {
	const Point3f& p0 = mesh->p[v[0]];
	const Point3f& p1 = mesh->p[v[1]];
	const Point3f& p2 = mesh->p[v[2]];

	// Translate vertices based on ray origin, so ray is originates from (0,0,0)
	Point3f p0t = p0 - Vector3f(ray.o);
	Point3f p1t = p1 - Vector3f(ray.o);
	Point3f p2t = p2 - Vector3f(ray.o);

	// Permute vertices and ray so largest ray dimension in z-axis (x and y axis arbitrary) --> +z direction
	int kz = maxDimension(abs(ray.d));
	int kx = kz + 1;
	if (kx == 3) {
		kx = 0;
	}
	int ky = kx + 1;
	if (ky == 3) {
		ky = 0;
	}

	Vector3f d = permute(ray.d, kx, ky, kz);
	p0t = permute(p0t, kx, ky, kz);
	p1t = permute(p1t, kx, ky, kz);
	p2t = permute(p2t, kx, ky, kz);

	// Shear to align ray along z-axis
	float Sx = -d.x / d.z;
	float Sy = -d.y / d.z;
	float Sz = 1.f / d.z;
	p0t.x += Sx * p0t.z;
	p0t.y += Sy * p0t.z;
	p1t.x += Sx * p1t.z;
	p1t.y += Sy * p1t.z;
	p2t.x += Sx * p2t.z;
	p2t.y += Sy * p2t.z;

	// Now with transformation, just have to check if point (0,0) is inside triangle in x,y-plane
	float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
	float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
	float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

	// Fall back to double precision test at triangle edges
	if (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f) {
		double p2txp1ty = (double)p2t.x * (double)p1t.y;
		double p2typ1tx = (double)p2t.y * (double)p1t.x;
		e0 = (float)(p2typ1tx - p2txp1ty);
		double p0txp2ty = (double)p0t.x * (double)p2t.y;
		double p0typ2tx = (double)p0t.y * (double)p2t.x;
		e1 = (float)(p0typ2tx - p0txp2ty);
		double p1txp0ty = (double)p1t.x * (double)p0t.y;
		double p1typ0tx = (double)p1t.y * (double)p0t.x;
		e2 = (float)(p1typ0tx - p1txp0ty);
	}

	// Triangle edge and determinant tests
	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0)) {
		return false;
	}
	float det = e0 + e1 + e2;
	if (det == 0) {
		return false;
	}

	// Compute scaled hit distance t and test
	p0t.z *= Sz;
	p1t.z *= Sz;
	p2t.z *= Sz;
	float t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
	if (det < 0 && (t_scaled >= 0 || t_scaled < ray.t_max * det)) {
		return false;
	} else if (det > 0 && (t_scaled <= 0 || t_scaled > ray.t_max * det)) {
		return false;
	}

	// Compute barycentric coords and t value
	float inv_det = 1 / det;
	float b0 = e0 * inv_det;
	float b1 = e1 * inv_det;
	float b2 = e2 * inv_det;
	float t = t_scaled * inv_det;

	// Ensure that computed triangle $t$ is conservatively greater than zero

	// Compute $\delta_z$ term for triangle $t$ error bounds
	float max_zt = maxComponent(abs(Vector3f(p0t.z, p1t.z, p2t.z)));
	float delta_z = gamma(3) * max_zt;

	// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
	float max_xt = maxComponent(abs(Vector3f(p0t.x, p1t.x, p2t.x)));
	float max_yt = maxComponent(abs(Vector3f(p0t.y, p1t.y, p2t.y)));
	float delta_x = gamma(5) * (max_xt + max_zt);
	float delta_y = gamma(5) * (max_yt + max_zt);

	// Compute $\delta_e$ term for triangle $t$ error bounds
	float deltaE =
			2 * (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);

	// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
	float max_e = maxComponent(abs(Vector3f(e0, e1, e2)));
	float delta_t = 3 *
			(gamma(3) * max_e * max_zt + deltaE * max_zt + delta_z * max_e) *
			std::abs(inv_det);
	if (t <= delta_t) {
		return false;
	}

	// Now the partial derivatives
	Vector3f dpdu, dpdv;
	Point2f uv[3];
	getUVs(uv);
	Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
	Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
	float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
	bool is_singular = std::abs(determinant) < 1e-8f;
	if (!is_singular) {
		float inv_determinant = 1 / determinant;
		dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * inv_determinant;
		dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * inv_determinant;
	}
	if (is_singular || cross(dpdu, dpdv).lengthSquared() == 0) {
		// if matrix is singular, choose arbitrary coord system
		Vector3f ng = cross(p2 - p0, p1 - p0);
		if (ng.lengthSquared() == 0) {
			return false;
		}
		coordinateSystem(normalize(ng), &dpdu, &dpdv);
	}

	// Compute error bounds for triangle intersection
	float x_abs_sum = (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
	float y_abs_sum = (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
	float z_abs_sum = (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
	Vector3f p_error = gamma(7) * Vector3f(x_abs_sum, y_abs_sum, z_abs_sum);

	Point3f p_hit = b0 * p0 + b1 * p1 + b2 * p2;
	Point2f uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

	// Test intersection against alpha mask, if valid
	if (test_alpha_texture && mesh->alpha_mask) {
		SurfaceInteraction isect_local(p_hit, Vector3f(0, 0, 0), uv_hit, Vector3f(0, 0, 0),
				dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);
		// TODO: implement texture class
//		if (mesh->alpha_mask->evaluate(isect_local) == 0) {
//			return false;
//		}
	}

	*isect = SurfaceInteraction(p_hit, p_error, uv_hit, -ray.d,
			dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);
	isect->n = isect->shading.n = Normal3f(normalize(cross(dp02, dp12))); ///< calculate surface normals manually using the edge vectors
	if (mesh->n || mesh->s) {
		// If shading normals provided, interpolate between them for smooth shading effect
		Normal3f ns;
		if (mesh->n) {
			ns = normalize(b0 * mesh->n[v[0]] + b1 * mesh->n[v[1]] + b2 * mesh->n[v[2]]);
			if (ns.lengthSquared() > 0) {
				ns = normalize(ns);
			} else {
				ns = isect->n;
			}
		} else {
			ns = isect->n;
		}

		Vector3f ss;
		if (mesh->s) {
			ss = normalize(b0 * mesh->s[v[0]] + b1 * mesh->s[v[1]] + b2 * mesh->s[v[2]]);
			if (ss.lengthSquared() > 0) {
				ss = normalize(ss);
			} else {
				ss = normalize(isect->dpdu);
			}
		} else {
			ss = normalize(isect->dpdu);
		}

		// Compute bitangent and adjust ss accordingly (if not perfectly orthogonal)
		Vector3f ts = cross(ss, ns);
		if (ts.lengthSquared() > 0) {
			ts = normalize(ts);
			ss = cross(ts, ns);
		} else {
			coordinateSystem(Vector3f(ns), &ss, &ts);
		}

		// Compute dndu and dndv for shading geometry
		Normal3f dndu, dndv;
		if (mesh->n) {
			// Compute deltas for triangle partial derivatives of normal
			Vector2f duv02 = uv[0] - uv[2];
			Vector2f duv12 = uv[1] - uv[2];
			Normal3f dn1 = mesh->n[v[0]] - mesh->n[v[2]];
			Normal3f dn2 = mesh->n[v[1]] - mesh->n[v[2]];
			float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
			bool is_singular = std::abs(determinant) < 1e-8f;
			if (is_singular) {
				// Can use same arbitrary coord system as for dpdu, dpdv to still get reasonable result
				Vector3f dn = cross(Vector3f(mesh->n[v[2]] - mesh->n[v[0]]),
						Vector3f(mesh->n[v[1]] - mesh->n[v[0]]));
				if (dn.lengthSquared() == 0) {
					dndu = dndv = Normal3f(0, 0, 0);
				} else {
					Vector3f dnu, dnv;
					coordinateSystem(dn, &dnu, &dnv);
					dndu = Normal3f(dnu);
					dndv = Normal3f(dnv);
				}
			} else {
				float inv_determinant = 1 / determinant;
				dndu = (duv12[1] * dn1 - duv02[1] * dn2) * inv_determinant;
				dndv = (-duv12[0] * dn1 + duv02[0] * dn2) * inv_determinant;
			}
		} else {
			dndu = dndv = Normal3f(0, 0, 0);
		}
		if (reverse_orientation) {
			ts = -ts;
		}
		isect->setShadingGeometry(ss, ts, dndu, dndv, true);
	}

	*t_hit = t;
	return true;
}

CIVET_CPU_GPU
bool Triangle::intersectP(const Ray& ray, bool test_alpha_texture) const {
	const Point3f& p0 = mesh->p[v[0]];
	const Point3f& p1 = mesh->p[v[1]];
	const Point3f& p2 = mesh->p[v[2]];

	// Translate vertices based on ray origin, so ray is originates from (0,0,0)
	Point3f p0t = p0 - Vector3f(ray.o);
	Point3f p1t = p1 - Vector3f(ray.o);
	Point3f p2t = p2 - Vector3f(ray.o);

	// Permute vertices and ray so largest ray dimension in z-axis (x and y axis arbitrary) --> +z direction
	int kz = maxDimension(abs(ray.d));
	int kx = kz + 1;
	if (kx == 3) {
		kx = 0;
	}
	int ky = kx + 1;
	if (ky == 3) {
		ky = 0;
	}

	Vector3f d = permute(ray.d, kx, ky, kz);
	p0t = permute(p0t, kx, ky, kz);
	p1t = permute(p1t, kx, ky, kz);
	p2t = permute(p2t, kx, ky, kz);

	// Shear to align ray along z-axis
	float Sx = -d.x / d.z;
	float Sy = -d.y / d.z;
	float Sz = 1.f / d.z;
	p0t.x += Sx * p0t.z;
	p0t.y += Sy * p0t.z;
	p1t.x += Sx * p1t.z;
	p1t.y += Sy * p1t.z;
	p2t.x += Sx * p2t.z;
	p2t.y += Sy * p2t.z;


	// Now with transformation, just have to check if point (0,0) is inside triangle in x,y-plane
	float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
	float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
	float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

	// Fall back to double precision test at triangle edges
	if (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f) {
		double p2txp1ty = (double)p2t.x * (double)p1t.y;
		double p2typ1tx = (double)p2t.y * (double)p1t.x;
		e0 = (float)(p2typ1tx - p2txp1ty);
		double p0txp2ty = (double)p0t.x * (double)p2t.y;
		double p0typ2tx = (double)p0t.y * (double)p2t.x;
		e1 = (float)(p0typ2tx - p0txp2ty);
		double p1txp0ty = (double)p1t.x * (double)p0t.y;
		double p1typ0tx = (double)p1t.y * (double)p0t.x;
		e2 = (float)(p1typ0tx - p1txp0ty);
	}

	// Triangle edge and determinant tests
	if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0)) {
		return false;
	}
	float det = e0 + e1 + e2;
	if (det == 0) {
		return false;
	}

	// Compute scaled hit distance t and test
	p0t.z *= Sz;
	p1t.z *= Sz;
	p2t.z *= Sz;
	float t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
	if (det < 0 && (t_scaled >= 0 || t_scaled < ray.t_max * det)) {
		return false;
	} else if (det > 0 && (t_scaled <= 0 || t_scaled > ray.t_max * det)) {
		return false;
	}

	// Compute barycentric coords and t value
	float inv_det = 1 / det;
	float b0 = e0 * inv_det;
	float b1 = e1 * inv_det;
	float b2 = e2 * inv_det;
	float t = t_scaled * inv_det;

	// Ensure that computed triangle $t$ is conservatively greater than zero

	// Compute $\delta_z$ term for triangle $t$ error bounds
	float max_zt = maxComponent(abs(Vector3f(p0t.z, p1t.z, p2t.z)));
	float delta_z = gamma(3) * max_zt;

	// Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
	float max_xt = maxComponent(abs(Vector3f(p0t.x, p1t.x, p2t.x)));
	float max_yt = maxComponent(abs(Vector3f(p0t.y, p1t.y, p2t.y)));
	float delta_x = gamma(5) * (max_xt + max_zt);
	float delta_y = gamma(5) * (max_yt + max_zt);

	// Compute $\delta_e$ term for triangle $t$ error bounds
	float deltaE =
			2 * (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);

	// Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
	float max_e = maxComponent(abs(Vector3f(e0, e1, e2)));
	float delta_t = 3 *
			(gamma(3) * max_e * max_zt + deltaE * max_zt + delta_z * max_e) *
			std::abs(inv_det);
	if (t <= delta_t) {
		return false;
	}

	// Now the partial derivatives
	Vector3f dpdu, dpdv;
	Point2f uv[3];
	getUVs(uv);
	Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
	Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
	float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
	bool is_singular = std::abs(determinant) < 1e-8f;
	if (!is_singular) {
		float inv_determinant = 1 / determinant;
		dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * inv_determinant;
		dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * inv_determinant;
	}
	if (is_singular || cross(dpdu, dpdv).lengthSquared() == 0) {
		// if matrix is singular, choose arbitrary coord system
		Vector3f ng = cross(p2 - p0, p1 - p0);
		if (ng.lengthSquared() == 0) {
			return false;
		}
		coordinateSystem(normalize(ng), &dpdu, &dpdv);
	}

	// Compute error bounds for triangle intersection
	float x_abs_sum = (std::abs(b0 * p0.x) + std::abs(b1 * p1.x) + std::abs(b2 * p2.x));
	float y_abs_sum = (std::abs(b0 * p0.y) + std::abs(b1 * p1.y) + std::abs(b2 * p2.y));
	float z_abs_sum = (std::abs(b0 * p0.z) + std::abs(b1 * p1.z) + std::abs(b2 * p2.z));
	Vector3f p_error = gamma(7) * Vector3f(x_abs_sum, y_abs_sum, z_abs_sum);

	Point3f p_hit = b0 * p0 + b1 * p1 + b2 * p2;
	Point2f uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

	// Test intersection against alpha mask, if valid
	if (test_alpha_texture && mesh->alpha_mask) {
		SurfaceInteraction isect_local(p_hit, Vector3f(0, 0, 0), uv_hit, Vector3f(0, 0, 0),
				dpdu, dpdv, Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this);
		// TODO: also implement texture class
//		if (mesh->alpha_mask->evaluate(isect_local) == 0) {
//			return false;
//		}
	}

	return true;
}

CIVET_CPU_GPU
float Triangle::area() const {
	const Point3f& p0 = mesh->p[v[0]];
	const Point3f& p1 = mesh->p[v[1]];
	const Point3f& p2 = mesh->p[v[2]];
	return 0.5 * cross(p1 - p0, p2 - p0).length();
}

} // namespace civet