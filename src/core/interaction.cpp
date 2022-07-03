#include "interaction.h"

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

	if (shape && (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
		n *= -1;
		shading.n *= -1;
	}
}

CIVET_CPU_GPU
void SurfaceInteraction::setShadingGeometry(const Vector3f& s_dpdu,const Vector3f& s_dpdv, const Normal3f& s_dndu, const Normal3f& s_dndv, bool orientation_is_authoritative) {
	shading.n = (Normal3f)normalize(cross(s_dpdu, s_dpdv));
	if (shape && (shape->reverseOrientation ^ shape->transformSwapsHandedness)) {
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

} // namespace civet