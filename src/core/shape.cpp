#include <core/shape.h>

namespace civet {

Shape::~Shape() {
}

Shape::Shape(const civet::Transform* otw, const civet::Transform* wto, bool _reverse_orientation) :
		object_to_world(otw), world_to_object(wto), reverse_orientation(_reverse_orientation), transform_swaps_handedness(otw->swapsHandedness()) {
}

Bounds3f Shape::worldBound() const {
	return (*object_to_world)(objectBound());
}

} // namespace civet
