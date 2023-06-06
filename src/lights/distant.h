#ifndef CIVET_LIGHTS_DISTANT_H
#define CIVET_LIGHTS_DISTANT_H

#include <core/civet.h>
#include <core/light.h>
#include <core/scene.h>
#include <core/shape.h>

namespace civet {

class DistantLight : public Light {
public:
	DistantLight(const Transform& ltw, const Spectrum& L, const Vector3f& w) :
			Light((int)LightFlags::DeltaDirection, ltw, MediumInterface()), L(L), w_light(normalize(w)) {}

	void preprocess(const Scene& scene) override {
		scene.worldBound().boundingSphere(&world_center, &world_radius);
	}

	Spectrum sample_Li(const Interaction &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTester *vis) const override;

	float pdf_Li(const Interaction &ref, const Vector3f &wi) const override {
		return 0;
	}

	Spectrum power() const override;

private:
	const Spectrum L;
	const Vector3f w_light;
	Point3f world_center;
	float world_radius;
};

} // namespace civet

#endif // CIVET_LIGHTS_DISTANT_H
