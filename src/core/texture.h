#ifndef CIVET_TEXTURE_H
#define CIVET_TEXTURE_H

#include <core/civet.h>
#include <core/geometry/transform.h>

namespace civet {

class TextureMapping2D {
public:
	virtual ~TextureMapping2D();
	virtual Point2f map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const = 0;
};

class UVMapping2D : public TextureMapping2D {
public:
	UVMapping2D(float su = 1, float sv = 1, float du = 0, float dv = 0);
	Point2f map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const override;

private:
	const float su, sv, du, dv;
};

class SphericalMapping2D : public TextureMapping2D {
public:
	SphericalMapping2D(const Transform& wtt) :
			world_to_texture(wtt) {}

	Point2f map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const override;

private:
	Point2f sphere(const Point3f& p) const;
	const Transform world_to_texture;
};

class PlanarMapping2D : public TextureMapping2D {
public:
	PlanarMapping2D(const Vector3f& vs, const Vector3f& vt, float ds = 0, float dt = 0) :
			vs(vs), vt(vt), ds(ds), dt(dt) {}

	Point2f map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const override;

private:
	const Vector3f vs, vt;
	const float ds, dt;
};

class TextureMapping3D {
public:
	virtual Point3f map(const SurfaceInteraction& si, Vector3f* dpdx, Vector3f* dpdy) const = 0;
};

class TransformMapping3D : public TextureMapping3D {
public:
	TransformMapping3D(const Transform& wtt) :
			world_to_texture(wtt) {}

	Point3f map(const SurfaceInteraction &si, Vector3f *dpdx, Vector3f *dpdy) const override;

private:
	const Transform world_to_texture;
};

template <typename T>
class Texture {
public:
	virtual T evaluate(const SurfaceInteraction&) const = 0;
	virtual ~Texture() {}
};

float lanczos(float, float tau = 2);

} // namespace civet

#endif // CIVET_TEXTURE_H
