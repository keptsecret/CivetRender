#include <core/interaction.h>

#include <core/shape.h>

namespace civet {

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

} // namespace civet