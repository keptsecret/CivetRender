#ifndef CIVET_TRIANGLE_H
#define CIVET_TRIANGLE_H

#include <core/shape.h>
#include <core/mesh.h>
#include <vector>

namespace civet {

class Triangle : public Shape {
public:
	CIVET_CPU_GPU
	Triangle(const Transform* otw, const Transform* wto, bool _reverse_orientation, const std::shared_ptr<TriangleMesh> _mesh, int tri_num) :
			Shape(otw, wto, _reverse_orientation), mesh(_mesh) {
		v = &(mesh->vertex_indices[3 * tri_num]);
	}

	CIVET_CPU_GPU
	Bounds3f objectBound() const override;
	CIVET_CPU_GPU
	Bounds3f worldBound() const override;

	CIVET_CPU_GPU
	bool intersect(const Ray &ray, float *t_hit, SurfaceInteraction *isect, bool test_alpha_texture = true) const override;
	CIVET_CPU_GPU
	bool intersectP(const Ray &ray, bool test_alpha_texture = true) const override;
	CIVET_CPU_GPU
	float area() const override;

private:
	CIVET_CPU_GPU
	void getUVs(Point2f uv[3]) const {
		if (mesh->uv) {
			uv[0] = mesh->uv[v[0]];
			uv[1] = mesh->uv[v[1]];
			uv[2] = mesh->uv[v[2]];
		} else {
			uv[0] = Point2f(0, 0);
			uv[1] = Point2f(1, 0);
			uv[2] = Point2f(1, 1);
		}
	}

	std::shared_ptr<TriangleMesh> mesh;
	const int* v;
};

} // namespace civet

#endif // CIVET_TRIANGLE_H
