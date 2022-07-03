#ifndef CIVET_INTERACTION_H
#define CIVET_INTERACTION_H

#include "civet.h"
#include "vecmath.h"

namespace civet {

struct Interaction {
	Interaction() :
			time(0) {}

	CIVET_CPU_GPU
	Interaction(const Point3f& _p, const Normal3f& _n, const Vector3f& _pe, const Vector3f& _wo, float t, const MediumInterface& mi) :
			p(_p), time(t), p_error(_pe), wo(_wo), n(_n), medium_interface(mi) {}

	CIVET_CPU_GPU
	bool isSurfaceInteraction() const { return n != Normal3f(); }

	Point3f p;
	float time;
	Vector3f p_error;
	Vector3f wo; /// negative ray direction, is (0,0,0) when outgoing direction doesn't apply
	Normal3f n;
	MediumInterface medium_interface;
};

class SurfaceInteraction : public Interaction {
public:
	CIVET_CPU_GPU
	SurfaceInteraction() {}

	CIVET_CPU_GPU
	SurfaceInteraction(const Point3f& _p, const Vector3f& _pe, const Point2f& _uv, const Vector3f& _wo,
			const Vector3f& _dpdu, const Vector3f& _dpdv, const Normal3f& _dndu, const Normal3f& _dndv,
			float t, const Shape* sh, int fi = 0);

	CIVET_CPU_GPU
	void setShadingGeometry(const Vector3f& s_dpdu,const Vector3f& s_dpdv, const Normal3f& s_dndu, const Normal3f& s_dndv, bool orientation_is_authoritative);

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
