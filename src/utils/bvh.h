#ifndef CIVET_BVH_H
#define CIVET_BVH_H

#include <core/civet.h>
#include <core/primitive.h>
#include <atomic>

namespace civet {

struct BVHPrimitiveInfo;

struct BVHBuildNode;
struct MortonPrimitive;
struct LBVHTreelet;
struct LinearBVHNode;

class BVH : public Aggregate {
public:
	enum class SplitMethod { SAH,
		HLBVH,
		Middle,
		EqualCounts };

	BVH(const std::vector<std::shared_ptr<Primitive>>& ps, int max_prims = 1, SplitMethod sm = SplitMethod::SAH);

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
