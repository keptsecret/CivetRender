#ifndef CIVET_PRIMITIVE_H
#define CIVET_PRIMITIVE_H

#include <core/geometry/transform.h>
#include <core/civet.h>
#include <core/material.h>
#include <core/shape.h>

namespace civet {

class Primitive {
public:
	virtual Bounds3f worldBound() const = 0;

	virtual bool intersect(const Ray& r, SurfaceInteraction*) const = 0;
	virtual bool intersectP(const Ray& r) const = 0;

	virtual const AreaLight* getAreaLight() const = 0;
	virtual const Material* getMaterial() const = 0;

	virtual void computeScatteringFunctions(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const = 0;
};

class GeometricPrimitive : public Primitive {
public:
	GeometricPrimitive(const std::shared_ptr<Shape> _shape,
			const std::shared_ptr<Material> _material,
			const std::shared_ptr<AreaLight> _area_light,
			const MediumInterface* mi);

	virtual Bounds3f worldBound() const override;
	virtual bool intersect(const Ray& r, SurfaceInteraction* isect) const override;
	virtual bool intersectP(const Ray& r) const override;

	const AreaLight* getAreaLight() const override;
	const Material* getMaterial() const override;
	void computeScatteringFunction(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const override;

private:
	std::shared_ptr<Shape> shape;
	std::shared_ptr<Material> material;
	std::shared_ptr<AreaLight> area_light;
	const MediumInterface* medium_interface; // TODO: turn back into instance when done
};

class TransformedPrimitve : public Primitive {
public:
	TransformedPrimitve(std::shared_ptr<Primitive> _primitive, const AnimatedTransform& ptw);

	virtual Bounds3f worldBound() const override {
		return primitive_to_world.motionBounds(primitive->worldBound());
	}
	virtual bool intersect(const Ray& r, SurfaceInteraction* isect) const override;
	virtual bool intersectP(const Ray& r) const override;

	const AreaLight* getAreaLight() const override { return nullptr; }
	const Material* getMaterial() const override { return nullptr; }
	void computeScatteringFunction(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const override {
		printf("ERROR::TransformedPrimitive: computeScatteringFunctions shouldn't be called");
	}

private:
	std::shared_ptr<Primitive> primitive;
	const AnimatedTransform primitive_to_world;
};

class Aggregate : public Primitive {
	const AreaLight* getAreaLight() const override {
		printf("ERROR::Aggregate: getAreaLight should've gone to a GeometricPrimitive");
		return nullptr;
	}
	const Material* getMaterial() const override {
		printf("ERROR::Aggregate: getMaterial should've gone to a GeometricPrimitive");
		return nullptr;
	}
	void computeScatteringFunction(SurfaceInteraction* isect, MemoryArena& arena, TransportMode mode, bool allow_multiple_lobes) const override {
		printf("ERROR::Aggregate: computeScatteringFunctions shouldn't be called");
	}
};

} // namespace civet

#endif // CIVET_PRIMITIVE_H
