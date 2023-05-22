#ifndef CIVET_LIGHT_POINT_H
#define CIVET_LIGHT_POINT_H

#include <core/light.h>
#include <core/spectrum.h>

namespace civet {

class PointLight : public Light {
public:
	PointLight(const Transform& ltw, const MediumInterface& mi, const Spectrum& I) :
			Light((int)LightFlags::DeltaPosition, ltw, mi), I(I) {}

	PointLight(const Transform& ltw, const MediumInterface& mi, Point3f p, const Spectrum& I) :
			Light((int)LightFlags::DeltaPosition, ltw, mi), pos_light(p), I(I) {}

	Spectrum sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTester *vis) const override;
	Spectrum power() const override;

private:
	const Point3f pos_light;
	const Spectrum I;
};

} // namespace civet

#endif // CIVET_LIGHT_POINT_H
