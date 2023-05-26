#include <core/texture.h>

#include <core/interaction.h>

namespace civet {

TextureMapping2D::~TextureMapping2D() {}

TextureMapping3D::~TextureMapping3D() {}

UVMapping2D::UVMapping2D(float su, float sv, float du, float dv) :
		su(su), sv(sv), du(du), dv(dv) {}

Point2f UVMapping2D::map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const {
	*dstdx = Vector2f(su * si.dudx, sv * si.dvdx);
	*dstdy = Vector2f(su * si.dudy, sv * si.dvdy);

	return Point2f(su * si.uv[0] + du, sv * si.uv[1] + dv);
}

Point2f SphericalMapping2D::map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const {
	Point2f st = sphere(si.p);

	const float delta = 0.1f;
	Point2f stdx = sphere(si.p + delta * si.dpdx);
	*dstdx = (stdx - st) / delta;
	Point2f stdy = sphere(si.p + delta * si.dpdy);
	*dstdy = (stdy - st) / delta;

	if ((*dstdx)[1] > 0.5f) {
		(*dstdx)[1] = 1 - (*dstdx)[1];
	} else if ((*dstdx)[1] < -0.5f) {
		(*dstdx)[1] = -((*dstdx)[1] + 1);
	}
	if ((*dstdy)[1] > 0.5f) {
		(*dstdy)[1] = 1 - (*dstdy)[1];
	} else if ((*dstdy)[1] < -.5f) {
		(*dstdy)[1] = -((*dstdy)[1] + 1);
	}

	return st;
}

Point2f SphericalMapping2D::sphere(const Point3f& p) const {
	Vector3f v = normalize(world_to_texture(p) - Point3f(0, 0, 0));
	float theta = sphericalTheta(v), phi = sphericalPhi(v);
	return Point2f(theta * InvPi, phi * InvPi);
}

Point2f PlanarMapping2D::map(const SurfaceInteraction& si, Vector2f* dstdx, Vector2f* dstdy) const {
	Vector3f v(si.p);
	*dstdx = Vector2f(dot(si.dpdx, vs), dot(si.dpdx, vt));
	*dstdy = Vector2f(dot(si.dpdy, vs), dot(si.dpdy, vt));
	return Point2f(ds + dot(v, vs), dt + dot(v, vt));
}

Point3f TransformMapping3D::map(const SurfaceInteraction& si, Vector3f* dpdx, Vector3f* dpdy) const {
	*dpdx = world_to_texture(si.dpdx);
	*dpdy = world_to_texture(si.dpdy);
	return world_to_texture(si.p);
}

float lanczos(float x, float tau) {
	x = std::abs(x);
	if (x < 1e-5f) {
		return 1;
	}
	if (x > 1.f) {
		return 0;
	}
	x *= Pi;
	float s = std::sin(x * tau) / (x * tau);
	float lanczos = std::sin(x) / x;
	return s * lanczos;
}

} // namespace civet