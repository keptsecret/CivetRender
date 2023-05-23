#ifndef CIVET_BVH_H
#define CIVET_BVH_H

#include <core/civet.h>
#include <core/primitive.h>
#include <atomic>

namespace civet {

struct BVHPrimitiveInfo {
	BVHPrimitiveInfo() {}
	BVHPrimitiveInfo(size_t prim_no, const Bounds3f& b) :
			primitive_num(prim_no), bounds(b), centroid(0.5f * b.p_min + 0.5f * b.p_max) {}

	size_t primitive_num;
	Bounds3f bounds;
	Point3f centroid;
};

struct BVHBuildNode {
	void initLeaf(int first, int n, const Bounds3f& b) {
		first_prim_offset = first;
		n_primitives = n;
		bounds = b;
		children[0] = children[1] = nullptr;
	}

	void initInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
		children[0] = c0;
		children[1] = c1;
		bounds = bUnion(c0->bounds, c1->bounds);
		split_axis = axis;
		n_primitives = 0;
	}

	Bounds3f bounds;
	BVHBuildNode* children[2];
	int split_axis, first_prim_offset, n_primitives;
};

struct MortonPrimitive {
	int primitive_index;
	uint32_t morton_code;
};

struct LBVHTreelet {
	int start_idx, n_primitives;
	BVHBuildNode* build_nodes;
};

struct LinearBVHNode {
	Bounds3f bounds;
	union {
		int primitives_offset; // for leaf node
		int second_child_offset; // for interior node
	};
	uint16_t n_primitives;
	uint8_t axis;
	uint8_t pad[1];
};

class BVH : public Aggregate {
public:
	enum class SplitMethod { SAH,
		HLBVH,
		Middle,
		EqualCounts };

	BVH(const std::vector<std::shared_ptr<Primitive>>& ps, int max_prims, SplitMethod sm);

	~BVH() { freeAligned(nodes); }

	Bounds3f worldBound() const override;
	bool intersect(const Ray& ray, SurfaceInteraction* isect) const;
	bool intersectP(const Ray& ray) const;

private:
	BVHBuildNode* recursiveBuild(MemoryArena& arena,
			std::vector<BVHPrimitiveInfo>& primitive_info,
			int start, int end, int* total_nodes,
			std::vector<std::shared_ptr<Primitive>>& ordered_prims);

	BVHBuildNode* HLBVHBuild(MemoryArena& arena,
			const std::vector<BVHPrimitiveInfo>& primitive_info,
			int* total_nodes, std::vector<std::shared_ptr<Primitive>>& ordered_prims);

	BVHBuildNode* emitLBVH(BVHBuildNode*& build_nodes, const std::vector<BVHPrimitiveInfo>& primitive_info,
			MortonPrimitive* morton_primitives, int n_primitives, int* total_nodes,
			std::vector<std::shared_ptr<Primitive>>& ordered_prims, std::atomic<int>* ordered_prims_offset, int bit_idx) const;

	BVHBuildNode* buildUpperSAH(MemoryArena& arena, std::vector<BVHBuildNode*>& treelet_roots, int start, int end, int* total_nodes) const;

	int flattenBVHTree(BVHBuildNode* node, int* offset);

	const int max_prims_in_node;
	const SplitMethod split_method;
	std::vector<std::shared_ptr<Primitive>> primitives;
	LinearBVHNode* nodes = nullptr;
};

} // namespace civet

#endif // CIVET_BVH_H
