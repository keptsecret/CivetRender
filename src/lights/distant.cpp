#include <lights/distant.h>

namespace civet {

Spectrum DistantLight::sample_Li(const Interaction& ref, const Point2f& u, Vector3f* wi, float* pdf, VisibilityTester* vis) const {
	*wi = w_light;
	*pdf = 1;
	Point3f p_outside = ref.p + w_light * (2 * world_radius);
	*vis = VisibilityTester(ref, Interaction(p_outside, ref.time, medium_interface));
	return L;
}

Spectrum DistantLight::power() const {
	return L * Pi * world_radius * world_radius;
}

} // namespace civet