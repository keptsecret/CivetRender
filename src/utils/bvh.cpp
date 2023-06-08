#include <utils/bvh.h>

#include <core/interaction.h>
#include <utils/parallel.h>
#include <algorithm>

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

// BVH utility functions
inline uint32_t leftShift3(uint32_t x) {
	if (x == (1 << 10)) {
		x--;
	}
	x = (x | (x << 16)) & 0b00000011000000000000000011111111;
	// x = ---- --98 ---- ---- ---- ---- 7654 3210
	x = (x | (x << 8)) & 0b00000011000000001111000000001111;
	// x = ---- --98 ---- ---- 7654 ---- ---- 3210
	x = (x | (x << 4)) & 0b00000011000011000011000011000011;
	// x = ---- --98 ---- 76-- --54 ---- 32-- --10
	x = (x | (x << 2)) & 0b00001001001001001001001001001001;
	// x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
	return x;
}

inline uint32_t encodeMorton3(const Vector3f& v) {
	return (leftShift3(v.z) << 2) | (leftShift3(v.y) << 1) | (leftShift3(v.x));
}

static void radixSort(std::vector<MortonPrimitive>* v) {
	std::vector<MortonPrimitive> temp_vec(v->size());
	CIVET_CONSTEXPR int bits_per_pass = 6;
	CIVET_CONSTEXPR int n_bits = 30;
	CIVET_CONSTEXPR int n_passes = n_bits / bits_per_pass;

	for (int pass = 0; pass < n_passes; pass++) {
		int low_bit = pass * bits_per_pass;
		std::vector<MortonPrimitive>& in = (pass & 1) ? temp_vec : *v;
		std::vector<MortonPrimitive>& out = (pass & 1) ? *v : temp_vec;

		CIVET_CONSTEXPR int n_buckets = 1 << bits_per_pass;
		int bucket_count[n_buckets] = { 0 };
		CIVET_CONSTEXPR int bit_mask = (1 << bits_per_pass) - 1;
		for (const auto& mp : in) {
			int bucket = (mp.morton_code >> low_bit) & bit_mask;
			++bucket_count[bucket];
		}

		int out_idx[n_buckets];
		out_idx[0] = 0;
		for (int i = 1; i < n_buckets; i++) {
			out_idx[i] = out_idx[i - 1] + bucket_count[i - 1];
		}

		for (const auto& mp : in) {
			int bucket = (mp.morton_code >> low_bit) & bit_mask;
			out[out_idx[bucket]++] = mp;
		}
	}

	if (n_passes & 1) {
		swapElem(*v, temp_vec);
	}
}

BVH::BVH(const std::vector<std::shared_ptr<Primitive>>& ps, int max_prims, BVH::SplitMethod sm) :
		max_prims_in_node(std::min(255, max_prims)), split_method(sm), primitives(std::move(ps)) {
	if (primitives.empty()) {
		return;
	}

	// Build BVH from primitives
	std::vector<BVHPrimitiveInfo> primitive_info(primitives.size());
	for (size_t i = 0; i < primitives.size(); i++) {
		primitive_info[i] = { i, primitives[i]->worldBound() };
	}

	MemoryArena arena(1024 * 1024);
	int total_nodes = 0;
	std::vector<std::shared_ptr<Primitive>> ordered_prims;
	ordered_prims.reserve(primitives.size());
	BVHBuildNode* root;
	if (split_method == SplitMethod::HLBVH) {
		root = HLBVHBuild(arena, primitive_info, &total_nodes, ordered_prims);
	} else {
		root = recursiveBuild(arena, primitive_info, 0, primitives.size(), &total_nodes, ordered_prims);
	}
	primitives.swap(ordered_prims);

	primitive_info.resize(0);
	printf("BVH created with %d nodes for %d primitives (%.2f MB), arena allocated %.2f MB\n",
			total_nodes, int(primitives.size()), float(total_nodes * sizeof(LinearBVHNode)) / (1024 * 1024),
			float(arena.totalAllocated()) / (1024 * 1024));
	nodes = allocAligned<LinearBVHNode>(total_nodes);
	int offset = 0;
	flattenBVHTree(root, &offset);
	if (offset != total_nodes) {
		printf("Error::BVH::BVH: missing nodes, should have %d, found %d\n", total_nodes, offset);
	}
}

struct BucketInfo {
	int count = 0;
	Bounds3f bounds;
};

BVHBuildNode* BVH::recursiveBuild(MemoryArena& arena,
		std::vector<BVHPrimitiveInfo>& primitive_info,
		int start, int end, int* total_nodes,
		std::vector<std::shared_ptr<Primitive>>& ordered_prims) {
	BVHBuildNode* node = arena.alloc<BVHBuildNode>();
	(*total_nodes)++;
	Bounds3f bounds;
	for (int i = start; i < end; i++) {
		bounds = bUnion(bounds, primitive_info[i].bounds);
	}

	int n_primitives = end - start;
	if (n_primitives == 1) {
		// create leaf node
		int first_prim_offset = ordered_prims.size();
		for (int i = start; i < end; i++) {
			int prim_num = primitive_info[i].primitive_num;
			ordered_prims.push_back(primitives[prim_num]);
		}
		node->initLeaf(first_prim_offset, n_primitives, bounds);
		return node;
	} else {
		// choose axis to split along
		Bounds3f centroid_bounds;
		for (int i = start; i < end; i++) {
			centroid_bounds = bUnion(centroid_bounds, primitive_info[i].centroid);
		}
		int dim = centroid_bounds.maximumAxis();

		// partition primitives and build children
		int mid = (start + end) / 2;
		if (centroid_bounds.p_max[dim] == centroid_bounds.p_min[dim]) {
			// build leaf node
			int first_prim_offset = ordered_prims.size();
			for (int i = start; i < end; i++) {
				int prim_num = primitive_info[i].primitive_num;
				ordered_prims.push_back(primitives[prim_num]);
			}
			node->initLeaf(first_prim_offset, n_primitives, bounds);
			return node;
		} else {
			// partition based on split method
			switch (split_method) {
				case SplitMethod::Middle: {
					// partition through node's midpoint
					float p_mid = (centroid_bounds.p_min[dim] + centroid_bounds.p_max[dim]) / 2;
					BVHPrimitiveInfo* mid_ptr = std::partition(&primitive_info[start], &primitive_info[end - 1] + 1,
							[dim, p_mid](const BVHPrimitiveInfo& pi) {
								return pi.centroid[dim] < p_mid;
							});
					mid = mid_ptr - &primitive_info[0];

					// if primitives have large overlapping bounds, can fail to partition
					// then break and fall through to EqualCounts
					if (mid != start && mid != end) {
						break;
					}
					[[fallthrough]];
				}
				case SplitMethod::EqualCounts: {
					mid = (start + end) / 2;
					std::nth_element(&primitive_info[start], &primitive_info[mid], &primitive_info[end - 1] + 1,
							[dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
								return a.centroid[dim] < b.centroid[dim];
							});
					break;
				}
				case SplitMethod::SAH:
				default: {
					if (n_primitives <= 2) {
						// partition into equal sized subsets
						// since further computation with SAH isn't worth it
						mid = (start + end) / 2;
						std::nth_element(&primitive_info[start], &primitive_info[mid], &primitive_info[end - 1] + 1,
								[dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
									return a.centroid[dim] < b.centroid[dim];
								});
					} else {
						CIVET_CONSTEXPR int n_buckets = 12;
						BucketInfo buckets[n_buckets];

						for (int i = start; i < end; i++) {
							int b = n_buckets * centroid_bounds.offset(primitive_info[i].centroid)[dim];
							if (b == n_buckets) {
								b = n_buckets - 1;
							}
							buckets[b].count++;
							buckets[b].bounds = bUnion(buckets[b].bounds, primitive_info[i].bounds);
						}

						// SAH costs for splitting up to each bucket
						float cost[n_buckets - 1];
						for (int i = 0; i < n_buckets - 1; i++) {
							Bounds3f b0, b1;
							int count0 = 0, count1 = 0;
							for (int j = 0; j <= i; j++) {
								b0 = bUnion(b0, buckets[j].bounds);
								count0 += buckets[i].count;
							}
							for (int j = i + 1; j <= n_buckets; j++) {
								b1 = bUnion(b1, buckets[j].bounds);
								count1 += buckets[i].count;
							}
							cost[i] = 1.f + (count0 * b0.surfaceArea() + count1 * b1.surfaceArea()) / bounds.surfaceArea();
						}

						// find bucket to split with minimal SAH cost
						float min_cost = cost[0];
						int min_split_bucket = 0;
						for (int i = 1; i < n_buckets - 1; i++) {
							if (cost[i] < min_cost) {
								min_cost = cost[i];
								min_split_bucket = i;
							}
						}

						// create leaf if cost is low enough
						// else split primitives according to bucket
						float leaf_cost = n_primitives;
						if (n_primitives > max_prims_in_node || min_cost < leaf_cost) {
							BVHPrimitiveInfo* p_mid = std::partition(&primitive_info[start], &primitive_info[end - 1] + 1,
									[=](const BVHPrimitiveInfo& pi) {
										int b = n_buckets * centroid_bounds.offset(pi.centroid)[dim];
										if (b == n_buckets) {
											b = n_buckets - 1;
										}
										return b <= min_split_bucket;
									});
							mid = p_mid - &primitive_info[0];
						} else {
							// build leaf node
							int first_prim_offset = ordered_prims.size();
							for (int i = start; i < end; i++) {
								int prim_num = primitive_info[i].primitive_num;
								ordered_prims.push_back(primitives[prim_num]);
							}
							node->initLeaf(first_prim_offset, n_primitives, bounds);
							return node;
						}
					}
					break;
				}
			}

			node->initInterior(dim, recursiveBuild(arena, primitive_info, start, mid, total_nodes, ordered_prims),
					recursiveBuild(arena, primitive_info, mid, end, total_nodes, ordered_prims));
		}
	}
	return node;
}

BVHBuildNode* BVH::HLBVHBuild(MemoryArena& arena, const std::vector<BVHPrimitiveInfo>& primitive_info,
		int* total_nodes, std::vector<std::shared_ptr<Primitive>>& ordered_prims) {
	Bounds3f bounds;
	for (auto& pi : primitive_info) {
		bounds = bUnion(bounds, pi.centroid);
	}

	std::vector<MortonPrimitive> morton_prims(primitive_info.size());
	parallelFor([&](int i) {
		// initialize morton_prims[i]
		CIVET_CONSTEXPR int morton_bits = 10;
		CIVET_CONSTEXPR int morton_scale = 1 << morton_bits;
		morton_prims[i].primitive_index = primitive_info[i].primitive_num;
		Vector3f centroid_offset = bounds.offset(primitive_info[i].centroid);
		morton_prims[i].morton_code = encodeMorton3(centroid_offset * morton_scale);
	},
			primitive_info.size(), 512);

	// sort morton indices
	radixSort(&morton_prims);

	// create treelets at bottom of BVH
	std::vector<LBVHTreelet> trees_to_build;
	for (int start = 0, end = 1; end <= morton_prims.size(); end++) {
		uint32_t mask = 0b00111111111111000000000000000000;
		if (end == morton_prims.size() || ((morton_prims[start].morton_code & mask) != (morton_prims[end].morton_code & mask))) {
			// add entry to treelet vectors for this treelet
			int n_primitives = end - start;
			int max_bvh_nodes = 2 * n_primitives;
			BVHBuildNode* nodes = arena.alloc<BVHBuildNode>(max_bvh_nodes, false);
			trees_to_build.push_back({ start, n_primitives, nodes });
			start = end;
		}
	}

	// create treelets in parallel
	std::atomic<int> atomic_total(0), ordered_prims_offset(0);
	ordered_prims.resize(primitives.size());
	parallelFor([&](int i) {
		// generate ith treelet
		int nodes_created = 0;
		const int first_bit_idx = 29 - 12;
		LBVHTreelet& tr = trees_to_build[i];
		tr.build_nodes = emitLBVH(tr.build_nodes, primitive_info, &morton_prims[tr.start_idx],
				tr.n_primitives, &nodes_created, ordered_prims, &ordered_prims_offset, first_bit_idx);
		atomic_total += nodes_created;
	},
			trees_to_build.size());
	*total_nodes = atomic_total;

	// create and return SAH BVH from LBVH treelets
	std::vector<BVHBuildNode*> finished_treelets;
	finished_treelets.reserve(trees_to_build.size());
	for (auto& treelet : trees_to_build) {
		finished_treelets.push_back(treelet.build_nodes);
	}
	return buildUpperSAH(arena, finished_treelets, 0, finished_treelets.size(), total_nodes);
}

BVHBuildNode* BVH::emitLBVH(BVHBuildNode*& build_nodes, const std::vector<BVHPrimitiveInfo>& primitive_info,
		MortonPrimitive* morton_primitives, int n_primitives, int* total_nodes,
		std::vector<std::shared_ptr<Primitive>>& ordered_prims, std::atomic<int>* ordered_prims_offset, int bit_idx) const {
	if (bit_idx == -1 || n_primitives < max_prims_in_node) {
		// create leaf node
		(*total_nodes)++;
		BVHBuildNode* node = build_nodes++;
		Bounds3f bounds;
		int first_prim_offset = ordered_prims_offset->fetch_add(n_primitives);
		for (int i = 0; i < n_primitives; i++) {
			int primitive_idx = morton_primitives[i].primitive_index;
			ordered_prims[first_prim_offset + i] = primitives[primitive_idx];
			bounds = bUnion(bounds, primitive_info[primitive_idx].bounds);
		}
		node->initLeaf(first_prim_offset, n_primitives, bounds);
		return node;
	} else {
		int mask = 1 << bit_idx;

		// advance to next subtree if no splits
		if ((morton_primitives[0].morton_code & mask) == (morton_primitives[n_primitives - 1].morton_code & mask)) {
			return emitLBVH(build_nodes, primitive_info, morton_primitives, n_primitives, total_nodes, ordered_prims, ordered_prims_offset, bit_idx - 1);
		}

		// find split point for this dimension
		int search_start = 0, search_end = n_primitives - 1;
		while (search_start + 1 != search_end) {
			int mid = (search_start + search_end) / 2;
			if ((morton_primitives[search_start].morton_code & mask) == (morton_primitives[mid].morton_code & mask)) {
				search_start = mid;
			} else {
				search_end = mid;
			}
		}
		int split_offset = search_end;

		(*total_nodes)++;
		BVHBuildNode* node = build_nodes++;
		BVHBuildNode* lbvh[2] = {
			emitLBVH(build_nodes, primitive_info, morton_primitives,
					split_offset, total_nodes, ordered_prims, ordered_prims_offset, bit_idx - 1),
			emitLBVH(build_nodes, primitive_info, &morton_primitives[split_offset],
					n_primitives - split_offset, total_nodes, ordered_prims, ordered_prims_offset, bit_idx - 1)
		};
		int axis = bit_idx % 3;
		node->initInterior(axis, lbvh[0], lbvh[1]);
		return node;
	}
}

BVHBuildNode* BVH::buildUpperSAH(MemoryArena& arena, std::vector<BVHBuildNode*>& treelet_roots, int start, int end, int* total_nodes) const {
	int nNodes = end - start;
	if (nNodes == 1) {
		return treelet_roots[start];
	}
	(*total_nodes)++;
	BVHBuildNode* node = arena.alloc<BVHBuildNode>();

	Bounds3f bounds;
	for (int i = start; i < end; ++i) {
		bounds = bUnion(bounds, treelet_roots[i]->bounds);
	}

	Bounds3f centroid_bounds;
	for (int i = start; i < end; ++i) {
		Point3f centroid =
				(treelet_roots[i]->bounds.p_min + treelet_roots[i]->bounds.p_max) *
				0.5f;
		centroid_bounds = bUnion(centroid_bounds, centroid);
	}
	int dim = centroid_bounds.maximumAxis();

	CIVET_CONSTEXPR int n_buckets = 12;
	struct BucketInfo {
		int count = 0;
		Bounds3f bounds;
	};
	BucketInfo buckets[n_buckets];

	for (int i = start; i < end; ++i) {
		float centroid = (treelet_roots[i]->bounds.p_min[dim] + treelet_roots[i]->bounds.p_max[dim]) * 0.5f;
		int b = n_buckets * ((centroid - centroid_bounds.p_min[dim]) / (centroid_bounds.p_max[dim] - centroid_bounds.p_min[dim]));
		if (b == n_buckets) {
			b = n_buckets - 1;
		}
		buckets[b].count++;
		buckets[b].bounds = bUnion(buckets[b].bounds, treelet_roots[i]->bounds);
	}

	// Compute costs for splitting after each bucket
	float cost[n_buckets - 1];
	for (int i = 0; i < n_buckets - 1; ++i) {
		Bounds3f b0, b1;
		int count0 = 0, count1 = 0;
		for (int j = 0; j <= i; ++j) {
			b0 = bUnion(b0, buckets[j].bounds);
			count0 += buckets[j].count;
		}
		for (int j = i + 1; j < n_buckets; ++j) {
			b1 = bUnion(b1, buckets[j].bounds);
			count1 += buckets[j].count;
		}
		cost[i] = .125f +
				(count0 * b0.surfaceArea() + count1 * b1.surfaceArea()) /
						bounds.surfaceArea();
	}

	float min_cost = cost[0];
	int min_split_bucket = 0;
	for (int i = 1; i < n_buckets - 1; ++i) {
		if (cost[i] < min_cost) {
			min_cost = cost[i];
			min_split_bucket = i;
		}
	}

	BVHBuildNode** p_mid = std::partition(
			&treelet_roots[start], &treelet_roots[end - 1] + 1,
			[=](const BVHBuildNode* node) {
				float centroid = (node->bounds.p_min[dim] + node->bounds.p_max[dim]) * 0.5f;
				int b = n_buckets * ((centroid - centroid_bounds.p_min[dim]) / (centroid_bounds.p_max[dim] - centroid_bounds.p_min[dim]));
				if (b == n_buckets) {
					b = n_buckets - 1;
				}

				return b <= min_split_bucket;
			});
	int mid = p_mid - &treelet_roots[0];

	node->initInterior(
			dim, this->buildUpperSAH(arena, treelet_roots, start, mid, total_nodes),
			this->buildUpperSAH(arena, treelet_roots, mid, end, total_nodes));
	return node;
}

int BVH::flattenBVHTree(BVHBuildNode* node, int* offset) {
	LinearBVHNode* linear_node = &nodes[*offset];
	linear_node->bounds = node->bounds;
	int curr_offset = (*offset)++;
	if (node->n_primitives > 0) {
		linear_node->primitives_offset = node->first_prim_offset;
		linear_node->n_primitives = node->n_primitives;
	} else {
		linear_node->axis = node->split_axis;
		linear_node->n_primitives = 0;
		flattenBVHTree(node->children[0], offset);
		linear_node->second_child_offset = flattenBVHTree(node->children[1], offset);
	}
	return curr_offset;
}

Bounds3f BVH::worldBound() const {
	return nodes ? nodes[0].bounds : Bounds3f();
}

bool BVH::intersect(const Ray& ray, SurfaceInteraction* isect) const {
	if (!nodes) {
		return false;
	}

	bool hit = false;
	Vector3f inv_dir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dir_is_neg[3] = { inv_dir.x < 0, inv_dir.y < 0, inv_dir.z < 0 };

	int to_visit_offset = 0, current_node_idx = 0;
	int nodes_to_visit[64];
	while (true) {
		const LinearBVHNode* node = &nodes[current_node_idx];
		float t0, t1;
		if (node->bounds.intersectP(ray.o, ray.d, ray.t_max, &t0, &t1)) {
			if (node->n_primitives > 0) {
				// intersect ray with primitives in leaf node
				for (int i = 0; i < node->n_primitives; i++) {
					if (primitives[node->primitives_offset + i]->intersect(ray, isect)) {
						hit = true;
					}
				}
				if (to_visit_offset == 0) {
					break;
				}
				current_node_idx = nodes_to_visit[--to_visit_offset];
			} else {
				if (dir_is_neg[node->axis]){
					nodes_to_visit[to_visit_offset++] = current_node_idx + 1;
					current_node_idx = node->second_child_offset;
				} else {
					nodes_to_visit[to_visit_offset++] = node->second_child_offset;
					current_node_idx = current_node_idx + 1;
				}
			}
		} else {
			if (to_visit_offset == 0) {
				break;
			}
			current_node_idx = nodes_to_visit[--to_visit_offset];
		}
	}

	return hit;
}

bool BVH::intersectP(const Ray& ray) const {
	if (!nodes) {
		return false;
	}

	Vector3f inv_dir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
	int dir_is_neg[3] = { inv_dir.x < 0, inv_dir.y < 0, inv_dir.z < 0 };

	int to_visit_offset = 0, current_node_idx = 0;
	int nodes_to_visit[64];
	while (true) {
		const LinearBVHNode* node = &nodes[current_node_idx];
		float t0, t1;
		if (node->bounds.intersectP(ray.o, ray.d, ray.t_max, &t0, &t1)) {	// TODO: check other intersectP that doesn't work
			if (node->n_primitives > 0) {
				// intersect ray with primitives in leaf node
				for (int i = 0; i < node->n_primitives; i++) {
					if (primitives[node->primitives_offset + i]->intersectP(ray)) {
						return true;
					}
				}
				if (to_visit_offset == 0) {
					break;
				}
				current_node_idx = nodes_to_visit[--to_visit_offset];
			} else {
				if (dir_is_neg[node->axis]){
					nodes_to_visit[to_visit_offset++] = current_node_idx + 1;
					current_node_idx = node->second_child_offset;
				} else {
					nodes_to_visit[to_visit_offset++] = node->second_child_offset;
					current_node_idx = current_node_idx + 1;
				}
			}
		} else {
			if (to_visit_offset == 0) {
				break;
			}
			current_node_idx = nodes_to_visit[--to_visit_offset];
		}
	}

	return false;
}

} // namespace civet