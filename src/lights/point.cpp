#include <lights/point.h>

namespace civet {

Spectrum PointLight::sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, float* pdf, VisibilityTester* vis) const {
	*wi = normalize(pos_light - ref.p);
	*pdf = 1.f;
	*vis = VisibilityTester(ref, Interaction(pos_light, ref.time, medium_interface));
	return I / distanceSquared(pos_light, ref.p);
}

Spectrum PointLight::power() const {
	return 4 * Pi * I;
}

} // namespace civet