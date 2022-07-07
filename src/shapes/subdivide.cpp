#include <shapes/subdivide.h>

#include <shapes/triangle.h>
#include <map>
#include <set>

namespace civet {

struct SDVertex;
struct SDFace;

// Macros for looping convenience
#define NEXT(i) (((i) + 1) % 3)
#define PREV(i) (((i) + 2) % 3)

struct SDVertex {
	SDVertex(const Point3f& _p = Point3f(0, 0, 0)) :
			p(_p) {}

	/**
	 * Returns the number of adjacent vertices
	 * @return number of adjacent vertices
	 */
	int valence();

	/**
	 * Returns the ring of adjacent points around this vertex.\n
	 * If vertex is a boundary, then the adjacent boundary vertices will be the first and last elements of the pointer array
	 * @param p pointer array to store list of adjacent vertices
	 */
	void oneRing(Point3f* p);

	Point3f p;
	SDFace* start_face = nullptr;
	SDVertex* child = nullptr;
	bool regular = false, boundary = false;
};

struct SDFace {
	SDFace() {
		for (int i = 0; i < 3; i++) {
			v[i] = nullptr;
			f[i] = nullptr;
		}
		for (int i = 0; i < 4; i++) {
			children[i] = nullptr;
		}
	}

	int vnum(SDVertex* vert) const {
		for (int i = 0; i < 3; i++) {
			if (v[i] == vert) {
				return i;
			}
		}
		printf("ERROR::SDFace: Basic logic error in SDFace::vnum()");
		return -1;
	}

	SDVertex* otherVert(SDVertex* v0, SDVertex* v1) {
		for (int i = 0; i < 3; i++) {
			if (v[i] != v0 && v[i] != v1) {
				return v[i];
			}
			printf("ERROR::SDFace: Basic logic error in SDFace::otherVert()");
			return nullptr;
		}
	}

	SDFace* next_face(SDVertex* vert) { return f[vnum(vert)]; }
	SDFace* prev_face(SDVertex* vert) { return f[PREV(vnum(vert))]; }
	SDVertex* next_vertex(SDVertex* vert) { return v[NEXT(vnum(vert))]; }
	SDVertex* prev_vertex(SDVertex* vert) { return v[PREV(vnum(vert))]; }

	SDVertex* v[3];
	SDFace* f[3]; ///< indexed by: ith neighbor face is on edge from v[i] to v[(i+1)%3]
	SDFace* children[4];
};

struct SDEdge {
	SDEdge(SDVertex* v0 = nullptr, SDVertex* v1 = nullptr) {
		v[0] = std::min(v0, v1);
		v[1] = std::max(v0, v1);
		f[0] = f[1] = nullptr;
		f0edge_num = -1;
	}

	bool operator<(const SDEdge& e2) const {
		if (v[0] == e2.v[0]) {
			return v[1] < e2.v[1];
		}
		return v[0] < e2.v[0];
	}

	SDVertex* v[2];
	SDFace* f[2];
	int f0edge_num;
};

/**
 * Returns a new weighted position for a regular vertex based on the adjacent vertices around it.\n
 * Follows formula: (1 - valence * beta) * vertex_position + SUM(beta * adj_vertex_position)
 * @param vert vertex to adjust weight
 * @param beta beta weight
 * @return new position of the vertex
 */
static Point3f weightOneRing(SDVertex* vert, float beta) {
	int valence = vert->valence();
	Point3f* p_ring = ALLOCA(Point3f, valence);
	vert->oneRing(p_ring);
	Point3f p = (1 - valence * beta) * vert->p;
	for (int i = 0; i < valence; i++) {
		p += beta * p_ring[i];
	}
	return p;
}

/**
 * Returns a new weighted position for a boundary vertex based on ONLY the adjacent boundary vertices around it.\n
 * Follows formula: (1 - valence * beta) * vertex_position + SUM(beta * adj_vertex_position)
 * @param vert vertex to adjust weight
 * @param beta beta weight
 * @return new position of the vertex
 */
static Point3f weightBoundary(SDVertex* vert, float beta) {
	int valence = vert->valence();
	Point3f* p_ring = ALLOCA(Point3f, valence);
	vert->oneRing(p_ring);
	Point3f p = (1 - 2 * beta) * vert->p;
	p += beta * p_ring[0];
	p += beta * p_ring[valence - 1];
	return p;
}

inline int SDVertex::valence() {
	SDFace* f = start_face;
	if (!boundary) {
		// for interior vertices
		int nf = 1;
		while ((f = f->next_face(this)) != start_face) {
			nf++;
		}
		return nf;
	} else {
		// for boundary vertices
		int nf = 1;
		while ((f = f->next_face(this)) != nullptr) {
			nf++;
		}
		f = start_face;
		while ((f = f->prev_face(this)) != nullptr) {
			nf++;
		}
		return nf + 1;
	}
}

inline void SDVertex::oneRing(Point3f* p) {
	if (!boundary) {
		SDFace* face = start_face;
		do {
			*p++ = face->next_vertex(this)->p;
			face = face->next_face(this);
		} while (face != start_face);
	} else {
		SDFace* face = start_face;
		SDFace* f2;
		while ((f2 = face->next_face(this)) != nullptr) {
			face = f2;
		}
		*p++ = face->next_vertex(this)->p;
		do {
			*p++ = face->prev_vertex(this)->p;
			face = face->prev_face(this);
		} while (face != nullptr);
	}
}

inline float beta(int valence) {
	if (valence == 3) {
		return 3.0f / 16.0f;
	} else {
		return 3.0f / (8.0f * valence);
	}
}

inline float loopGamma(int valence) {
	return 1.0f / (valence + 3.0f / (8.0f * beta(valence)));
}

static std::vector<std::shared_ptr<Shape>> subdivide(const Transform* otw, const Transform* wto, bool _reverse_orientation,
		int n_levels, int n_indices, const int* v_idx, int n_verts, const Point3f* p) {
	std::vector<SDVertex*> vertices;
	std::vector<SDFace*> faces;

	std::unique_ptr<SDVertex[]> verts(new SDVertex[n_verts]);
	for (int i = 0; i < n_verts; i++) {
		verts[i] = SDVertex(p[i]);
		vertices.push_back(&verts[i]);
	}

	int n_faces = n_indices / 3;
	std::unique_ptr<SDFace[]> fs(new SDFace[n_faces]);
	for (int i = 0; i < n_faces; i++) {
		faces.push_back(&fs[i]);
	}

	// Set face to vertex pointers
	const int* vp = v_idx;
	for (int i = 0; i < n_faces; i++, vp += 3) {
		SDFace* f = faces[i];
		for (int j = 0; j < 3; j++) {
			SDVertex* v = vertices[vp[j]];
			f->v[j] = v;
			v->start_face = f;
		}
	}

	// Set neighbor face pointers in faces
	std::set<SDEdge> edges;
	for (int i = 0; i < n_faces; i++) {
		SDFace* f = faces[i];
		for (int edge_num = 0; edge_num < 3; edge_num++) {
			int v0 = edge_num, v1 = NEXT(edge_num);
			SDEdge e(f->v[v0], f->v[v1]);
			if (edges.find(e) == edges.end()) {
				// new edge
				e.f[0] = f;
				e.f0edge_num = edge_num;
				edges.insert(e);
			} else {
				// seen edge before
				e = *edges.find(e);
				e.f[0]->f[e.f0edge_num] = f;
				f->f[edge_num] = e.f[0];
				edges.erase(e);
			}
		}
	}

	// Determine if vertex is regular and if it is boundary
	// (see if we can follow adjacent faces and end up at same face)
	for (int i = 0; i < n_verts; i++) {
		SDVertex* v = vertices[i];
		SDFace* f = v->start_face;
		do {
			f = f->next_face(v);
		} while (f && f != v->start_face);
		v->boundary = (f == nullptr);
		if (!v->boundary && v->valence() == 6) {
			v->regular = true;
		} else if (v->boundary && v->valence() == 4) {
			v->regular = true;
		} else {
			v->regular = false;
		}
	}

	// Refine subdivision mesh into triangles
	std::vector<SDFace*> f = faces;
	std::vector<SDVertex*> v = vertices;
	MemoryArena arena;
	for (int i = 0; i < n_levels; i++) {
		std::vector<SDFace*> new_faces;
		std::vector<SDVertex*> new_vertices;

		for (auto vertex : v) {
			vertex->child = arena.alloc<SDVertex>();
			vertex->child->regular = vertex->regular;
			vertex->child->boundary = vertex->boundary;
			new_vertices.push_back(vertex->child);
		}

		for (auto face : f) {
			for (int k = 0; k < 4; k++) {
				face->children[k] = arena.alloc<SDFace>();
				new_faces.push_back(face->children[k]);
			}
		}

		// Update vertex positions
		for (auto vertex : v) {
			if (!vertex->boundary) {
				// apply one-ring rule for even vertex
				if (vertex->regular) {
					vertex->child->p = weightOneRing(vertex, 1.0f / 16.0f);
				} else {
					vertex->child->p = weightOneRing(vertex, beta(vertex->valence()));
				}
			} else {
				vertex->child->p = weightBoundary(vertex, 1.0f / 8.0f);
			}
		}

		// Compute new vertex positions
		std::map<SDEdge, SDVertex*> edge_verts;
		for (auto face : f) {
			for (int k = 0; k < 3; k++) {
				SDEdge edge(face->v[k], face->v[NEXT(k)]);
				SDVertex* vert = edge_verts[edge];
				if (!vert) {
					// initialize new vertex
					vert = arena.alloc<SDVertex>();
					new_vertices.push_back(vert);
					vert->regular = true;
					vert->boundary = (face->f[k] == nullptr);
					vert->start_face = face->children[3]; ///< guaranteed to be adjacent to center face

					// apply edge rules and compute position
					if (vert->boundary) {
						vert->p = 0.5 * edge.v[0]->p;
						vert->p += 0.5 * edge.v[1]->p;
					} else {
						vert->p = 3.0f / 8.0f * edge.v[0]->p;
						vert->p += 3.0f / 8.0f * edge.v[1]->p;
						vert->p += 1.0f / 8.0f * face->otherVert(edge.v[0], edge.v[1])->p;
						vert->p += 1.0f / 8.0f * face->f[k]->otherVert(edge.v[0], edge.v[1])->p;
					}
				}
			}
		}

		// Update mesh topology

		// update even vertex face pointers
		for (auto vertex : v) {
			int vert_num = vertex->start_face->vnum(vertex);
			vertex->child->start_face = vertex->start_face->children[vert_num];
		}
		// update face neighbor pointers
		for (auto face : f) {
			for (int j = 0; j < 3; j++) {
				// update pointers to neighbors of same parent
				face->children[3]->f[j] = face->children[NEXT(j)];
				face->children[j]->f[NEXT(j)] = face->children[3];

				// update pointers for neighbors of other parent face
				SDFace* f2 = face->f[j];
				face->children[j]->f[j] = f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
				f2 = face->f[PREV(j)];
				face->children[j]->f[PREV(j)] = f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
			}
		}
		// update face vertex pointers
		for (auto face : f) {
			for (int j = 0; j < 3; j++) {
				face->children[j]->v[j] = face->v[j]->child;

				SDVertex* vert = edge_verts[SDEdge(face->v[j], face->v[NEXT(j)])];
				face->children[j]->v[NEXT(j)] = vert;
				face->children[NEXT(j)]->v[j] = vert;
				face->children[3]->v[j] = vert;
			}
		}

		f = new_faces;
		v = new_vertices;
	}

	// Push vertices to limit surface
	std::unique_ptr<Point3f[]> p_limit(new Point3f[v.size()]);
	for (size_t i = 0; i < v.size(); i++) {
		if (v[i]->boundary) {
			p_limit[i] = weightBoundary(v[i], 1.0f / 5.0f);
		} else {
			p_limit[i] = weightOneRing(v[i], loopGamma(v[i]->valence()));
		}
	}
	for (size_t i = 0; i < v.size(); i++) {
		v[i]->p = p_limit[i];
	}

	// Compute vertex tangents on limit surface
	std::vector<Normal3f> Ns;
	Ns.reserve(v.size());
	std::vector<Point3f> p_ring(16, Point3f());
	for (auto vertex : v) {
		Vector3f S(0, 0, 0), T(0, 0, 0);
		int valence = vertex->valence();
		if (valence > (int)p_ring.size()) {
			p_ring.resize(valence);
		}
		vertex->oneRing(&p_ring[0]);
		if (!vertex->boundary) {
			for (int j = 0; j < valence; j++) {
				S += std::cos(2 * Pi * j / valence) * Vector3f(p_ring[j]);
				T += std::sin(2 * Pi * j / valence) * Vector3f(p_ring[j]);
			}
		} else {
			S = p_ring[valence - 1] - p_ring[0];
			if (valence == 2) {
				T = Vector3f(p_ring[0] + p_ring[1] - 2 * vertex->p);
			} else if (valence == 3) {
				T = p_ring[1] - vertex->p;
			} else if (valence == 4) { // regular
				T = Vector3f(-1 * p_ring[0] + 2 * p_ring[1] + 2 * p_ring[2] + -1 * p_ring[3] + -2 * vertex->p);
			} else {
				float theta = Pi / float(valence - 1);
				T = Vector3f(std::sin(theta) * (p_ring[0] + p_ring[valence - 1]));
				for (int k = 1; k < valence - 1; ++k) {
					float wt = (2 * std::cos(theta) - 2) * std::sin((k)*theta);
					T += Vector3f(wt * p_ring[k]);
				}
				T = -T;
			}
		}
		Ns.push_back(Normal3f(cross(S, T)));
	}

	// Create triangle mesh from subdivision mesh
	{
		size_t ntris = f.size();
		std::unique_ptr<int[]> verts(new int[3 * ntris]);
		int* vp = verts.get();
		size_t total_verts = v.size();
		std::map<SDVertex*, int> used_verts;
		for (size_t i = 0; i < total_verts; ++i) {
			used_verts[v[i]] = i;
		}
		for (size_t i = 0; i < ntris; ++i) {
			for (int j = 0; j < 3; ++j) {
				*vp = used_verts[f[i]->v[j]];
				++vp;
			}
		}
		return createTriangleMesh(otw, wto,
				_reverse_orientation, ntris, verts.get(),
				total_verts, p_limit.get(), nullptr, &Ns[0],
				nullptr, nullptr);
	}
}

} // namespace civet
