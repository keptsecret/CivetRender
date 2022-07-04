#ifndef CIVET_INTERACTION_H
#define CIVET_INTERACTION_H

#include <core/civet.h>
#include <core/vecmath.h>
#include <core/transform.h>
#include <core/shape.h>

namespace civet {

// TODO: remove when MediumInterface is implemented
class MediumInterface {

};

struct Interaction {
	Interaction() :
			time(0) {}

	CIVET_CPU_GPU
	Interaction(const Point3f& _p, const Normal3f& _n, const Vector3f& _pe, const Vector3f& _wo, float t, const MediumInterface* mi) :
			p(_p), time(t), p_error(_pe), wo(normalize(_wo)), n(_n), medium_interface(mi) {}

	CIVET_CPU_GPU
	Interaction(const Point3f& p, const Vector3f& wo, float time, const MediumInterface* mi) :
			p(p), time(time), wo(wo), medium_interface(mi) {}

	CIVET_CPU_GPU
	Interaction(const Point3f& p, float time, const MediumInterface* mi) :
			p(p), time(time), medium_interface(mi) {}

	CIVET_CPU_GPU
	bool isSurfaceInteraction() const { return n != Normal3f(); }

	Point3f p;
	float time;
	Vector3f p_error;
	Vector3f wo; /// negative ray direction, is (0,0,0) when outgoing direction doesn't apply
	Normal3f n;
	const MediumInterface* medium_interface = nullptr;
};

class SurfaceInteraction : public Interaction {
public:
	CIVET_CPU_GPU
	SurfaceInteraction() {}

	CIVET_CPU_GPU
	SurfaceInteraction(const Point3f& _p, const Vector3f& _pe, const Point2f& _uv, const Vector3f& _wo,
			const Vector3f& _dpdu, const Vector3f& _dpdv, const Normal3f& _dndu, const Normal3f& _dndv,
			float t, const Shape* sh, int fi = 0) :
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
	void setShadingGeometry(const Vector3f& s_dpdu, const Vector3f& s_dpdv, const Normal3f& s_dndu, const Normal3f& s_dndv, bool orientation_is_authoritative);

	Point2f uv;
	Vector3f dpdu, dpdv;
	Normal3f dndu, dndv;
	const Shape* shape = nullptr;

	// second instance for possibly perturbed values
	struct {
		Normal3f n;
		Vector3f dpdu, dpdv;
		Normal3f dndu, dndv;
	} shading;

	int face_index = 0;
};

} // namespace civet

#endif // CIVET_INTERACTION_H
