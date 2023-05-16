#ifndef CIVET_NODE_H
#define CIVET_NODE_H

#include <core/civet.h>
#include <core/geometry/transform.h>
#include <core/geometry/vecmath.h>

namespace civet {

enum NodeType {
	Model,
	Mesh,
	DirectionalLight,
	PointLight,
	SkyBox
};

struct TransformData {
	// local transform components
	Vector3f translation{ 0, 0, 0 };
	Vector3f rotation_vec{ 0, 0, 0 };
	Vector3f scale_vec{ 1, 1, 1 };

	TransformData() {
		updateTransform();
	}

	void updateTransform() {
		transform = translate(translation) * rotateX(rotation_vec.x) * rotateY(rotation_vec.y) * rotateZ(rotation_vec.z)
				* scale(scale_vec.x, scale_vec.y, scale_vec.z);
	}

	Transform transform;	// local transform
};

class Node {
public:
	Node(const std::string& n, const NodeType t, bool te = true) {
		name = n;
		type = t;
		transformEnabled = te;
	}

	virtual Transform getWorldTransform() {
		if (parent != nullptr) {
			return parent->getWorldTransform() * transform_data.transform;
		}

		return transform_data.transform;
	}

	virtual void updateWorldBounds() {
		bounds = Bounds3f();
		for (auto&& child : children) {
			child->updateWorldBounds();
			bounds = bUnion(bounds, child->getWorldBounds());
		}
	}

	virtual Bounds3f getWorldBounds() {
		return bounds;
	}

	bool isTransformEnabled() { return transformEnabled; }

	std::string name;
	NodeType type;

	Bounds3f bounds;	// global bounds
	TransformData transform_data;
	Node* parent = nullptr;
	std::vector<std::shared_ptr<Node>> children;

private:
	bool transformEnabled;
};

} // namespace civet

#endif // CIVET_NODE_H
