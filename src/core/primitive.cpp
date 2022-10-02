#include <core/primitive.h>

#include <core/interaction.h>

namespace civet {

GeometricPrimitive::GeometricPrimitive(const std::shared_ptr<Shape> _shape,
		const std::shared_ptr<Material> _material,
		const std::shared_ptr<AreaLight> _area_light,
		const MediumInterface* mi) :
		shape(_shape), material(_material), area_light(_area_light), medium_interface(mi) {}

Bounds3f GeometricPrimitive::worldBound() const {
	return shape->worldBound();
}

bool GeometricPrimitive::intersect(const Ray& r, SurfaceInteraction* isect) const {
	float t_hit;
	if (!shape->intersect(r, &t_hit, isect)) {
		return false;
	}
	r.t_max = t_hit;
	isect->primitive = (Primitive*)(this); ///< check if cast causes problems
	// TODO: initialize MediumInterface when it becomes implemented
	return true;
}

bool GeometricPrimitive::intersectP(const Ray& r) const {
	return shape->intersectP(r);
}

const AreaLight* GeometricPrimitive::getAreaLight() const {
	return area_light.get();
}

const Material* GeometricPrimitive::getMaterial() const {
	return material.get();
}

void GeometricPrimitive::computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const {
	if (material) {
		material->computeScatteringFunctions(isect, arena, mode, allow_multiple_lobes);
	}
}

TransformedPrimitve::TransformedPrimitve(std::shared_ptr<Primitive> _primitive, const AnimatedTransform& ptw) :
		primitive(_primitive), primitive_to_world(ptw) {}

bool TransformedPrimitve::intersect(const Ray& r, SurfaceInteraction* isect) const {
	Transform interpolated_ptw;
	primitive_to_world.interpolate(r.time, &interpolated_ptw);
	Ray ray = inverse(interpolated_ptw)(r);
	if (!primitive->intersect(ray, isect)) {
		return false;
	}
	r.t_max = ray.t_max;
	if (!interpolated_ptw.isIdentity()) {
		*isect = interpolated_ptw(*isect);
	}
	return true;
}

bool TransformedPrimitve::intersectP(const Ray& r) const {
	Transform interpolated_ptw;
	primitive_to_world.interpolate(r.time, &interpolated_ptw);
	Ray ray = inverse(interpolated_ptw)(r);
	return !primitive->intersectP(ray);
}

} // namespace civet