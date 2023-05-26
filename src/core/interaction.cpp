#include <core/interaction.h>

#include <core/geometry/transform.h>
#include <core/shape.h>
#include <core/primitive.h>

namespace civet {
CIVET_CPU_GPU
SurfaceInteraction::SurfaceInteraction(const Point3f& _p, const Vector3f& _pe, const Point2f& _uv, const Vector3f& _wo,
		const Vector3f& _dpdu, const Vector3f& _dpdv, const Normal3f& _dndu, const Normal3f& _dndv,
		float t, const Shape* sh, int fi) :
		Interaction(_p, Normal3f(normalize(cross(_dpdu, _dpdv))), _pe, _wo, t, nullptr),
		uv(_uv),
		dpdu(_dpdu),
		dpdv(_dpdv),
		dndu(_dndu),
		dndv(_dndv),
		shape(sh),
		face_index(fi) {
	// Initialize shading geometry from true geometry
	shading.n = n;
	shading.dpdu = dpdu;
	shading.dpdv = dpdv;
	shading.dndu = dndu;
	shading.dndv = dndv;

	if (shape && (shape->reverse_orientation ^ shape->transform_swaps_handedness)) {
		n *= -1;
		shading.n *= -1;
	}
}

CIVET_CPU_GPU
void SurfaceInteraction::setShadingGeometry(const Vector3f& s_dpdu, const Vector3f& s_dpdv, const Normal3f& s_dndu, const Normal3f& s_dndv, bool orientation_is_authoritative) {
	shading.n = (Normal3f)normalize(cross(s_dpdu, s_dpdv));
	if (shape && (shape->reverse_orientation ^ shape->transform_swaps_handedness)) {
		shading.n = -shading.n;
	}
	if (orientation_is_authoritative) {
		n = faceforward(n, shading.n);
	} else {
		n = faceforward(shading.n, n);
	}

	shading.dpdu = s_dpdu;
	shading.dpdv = s_dpdv;
	shading.dndu = s_dndu;
	shading.dndv = s_dndv;
}

CIVET_CPU_GPU
void SurfaceInteraction::computeDifferentials(const RayDifferential& ray) {
	if (ray.has_differentials) {
		float d = -dot(n, Vector3f(p.x, p.y, p.z));
		float tx = (-dot(n, Vector3f(ray.rx_origin)) - d) / dot(n, ray.rx_direction);
		Point3f px = ray.rx_origin + tx * ray.rx_direction;
		float ty = (-dot(n, Vector3f(ray.ry_origin)) - d) / dot(n, ray.ry_direction);
		Point3f py = ray.ry_origin + ty * ray.ry_direction;

		dpdx = px - p;
		dpdy = py - p;

		int dim[2];
		if (std::abs(n.x) > std::abs(n.y) && std::abs(n.x) > std::abs(n.z)) {
			dim[0] = 1, dim[1] = 2;
		} else if (std::abs(n.y) > std::abs(n.z)) {
			dim[0] = 0, dim[1] = 2;
		} else {
			dim[0] = 0, dim[1] = 1;
		}

		float A[2][2] = {{dpdu[dim[0]], dpdv[dim[0]]},
			{dpdu[dim[1]], dpdv[dim[1]]}};
		float Bx[2] = {px[dim[0]] - p[dim[0]], px[dim[1]] - p[dim[1]]};
		float By[2] = {py[dim[0]] - p[dim[0]], py[dim[1]] - p[dim[1]]};
		if (solveLinearSystem2x2(A, Bx, &dudx, &dvdx)) {
			dudx = dvdx = 0;
		}
		if (solveLinearSystem2x2(A, By, &dudy, &dvdy)) {
			dudy = dvdy = 0;
		}
	} else {
		dudx = dvdx = 0;
		dudy = dvdy = 0;
		dpdx = dpdy = Vector3f(0,0,0);
	}
}

CIVET_CPU_GPU
void SurfaceInteraction::computeScatteringFunctions(const RayDifferential& ray, MemoryArena& arena, bool allow_multiple_lobes, TransportMode mode) {
	computeDifferentials(ray);
	primitive->computeScatteringFunctions(this, arena, mode, allow_multiple_lobes);
}

} // namespace civet