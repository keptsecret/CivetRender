#ifndef CIVET_NODE_H
#define CIVET_NODE_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <core/geometry/transform.h>

namespace civet {

enum NodeType {
	Model,
	Mesh,
	DirectionalLight,
	PointLight
};

struct TransformData {
	// TODO: scene graph type transform data
	Vector3f position;
	Vector3f rotation;
	Vector3f scale;

	Transform transform;
};

class Node {
public:
	Node(const std::string& n, const NodeType t) {
		name = n;
		type = t;
	}

	virtual void updateBounds() {}
	virtual void setTransform(Transform t) {
		transformData.transform = t;
		updateBounds();
	}

	virtual Transform getWorldTransform() {
		if (parent != nullptr) {
			return parent->getWorldTransform() * transformData.transform;
		}

		return transformData.transform;
	}

	std::string name;
	NodeType type;

	Bounds3f bounds;
	TransformData transformData;
	Node* parent = nullptr;
};

} // namespace civet

#endif // CIVET_NODE_H
