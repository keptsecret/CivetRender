#ifndef CIVET_SHAPE_H
#define CIVET_SHAPE_H

#include <core/civet.h>
#include <core/transform.h>

namespace civet {

class Shape {
public:
	Shape(const Transform* otw, const Transform* wto, bool _reverse_orientation);

	virtual ~Shape();

	/**
	 * Returns a bounding box in the object space
	 * @return a Bounds3f object
	 */
	virtual Bounds3f objectBound() const = 0;

	/**
	 * Returns a bounding box in the world space, using the object_to_world transform by default
	 * @return a Bounds3f object
	 */
	virtual Bounds3f worldBound() const;

	/**
	 * Returns whether a ray intersection occurs with the shape and stores relevant information in provided structures
	 * @param ray incoming ray to test for intersection
	 * @param t_hit t value where the ray intersects the shape
	 * @param isect information on surface material and texture
	 * @param test_alpha_texture
	 * @return true if ray intersects shape
	 */
	virtual bool intersect(const Ray& ray, float* t_hit, SurfaceInteraction* isect, bool test_alpha_texture = true) const = 0;

	/**
	 * Returns whether a ray intersection occurs with the shape (acts only as a test, involves wasteful calculations)
	 * @param ray incoming ray to test for intersection
	 * @param test_alpha_texture
	 * @return true if ray intersects shape
	 */
	virtual bool intersectP(const Ray& ray, bool test_alpha_texture = true) const {
		return intersect(ray, nullptr, nullptr, test_alpha_texture);
	}

	virtual float area() const = 0;

	const Transform* object_to_world;
	const Transform* world_to_object;
	const bool reverse_orientation;
	const bool transform_swaps_handedness;
};

} // namespace civet

#endif // CIVET_SHAPE_H
