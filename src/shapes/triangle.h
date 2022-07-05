#ifndef CIVET_TRIANGLE_H
#define CIVET_TRIANGLE_H

#include <core/shape.h>
#include <vector>

namespace civet {

struct TriangleMesh {
	/**
	 * Initializes a triangle mesh
	 * @param otw object-to-world transformation for the mesh
	 * @param n_tris total number of triangles in the mesh
	 * @param v_idx pointer to array of vertex indices, where i-th triangle has vertices at elements (3i, 3i+1, 3i+2)
	 * @param n_verts total number of vertices in the mesh
	 * @param _p array of n_verts vertex positions
	 * @param _s (optional) array of tangent vectors, one per vertex for computing shading tangents
	 * @param _n (optional) array of normals, one per vertex, interpolated across triangle faces as shading normals
	 * @param _uv (optional) array of parametric u,v-values, one per vertex
	 * @param am (optional) alpha mask texture
	 */
	CIVET_CPU_GPU
	TriangleMesh(const Transform& otw, int n_tris, const int* v_idx, int n_verts,
			const Point3f* _p, const Vector3f* _s, const Normal3f* _n, const Point2f* _uv,
			std::shared_ptr<Texture<float>> am);

	const int n_triangles;
	const int n_vertices;
	std::vector<int> vertex_indices;
	std::unique_ptr<Point3f[]> p;
	std::unique_ptr<Normal3f[]> n;
	std::unique_ptr<Vector3f[]> s;
	std::unique_ptr<Point2f[]> uv;
	std::shared_ptr<Texture<float>> alpha_mask;
};

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

CIVET_CPU_GPU
std::vector<std::shared_ptr<Shape>> createTriangleMesh(const Transform* otw, const Transform* wto, bool _reverse_orientation,
		int n_tris, const int* v_idx, int n_verts,
		const Point3f* _p, const Vector3f* _s, const Normal3f* _n, const Point2f* _uv,
		std::shared_ptr<Texture<float>> am);

} // namespace civet

#endif // CIVET_TRIANGLE_H
