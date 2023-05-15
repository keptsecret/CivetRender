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
	Vector3f translation{ 0, 0, 0 };
	Vector3f rotation_vec{ 0, 0, 0 };
	Vector3f scale_vec{ 1, 1, 1 };

	TransformData() {
		transform = translate(translation) * rotateX(rotation_vec.x) * rotateY(rotation_vec.y) * rotateZ(rotation_vec.z) * scale(scale_vec.x, scale_vec.y, scale_vec.z);
	}

	void updateTransform() {
		transform = translate(translation) * rotateX(rotation_vec.x) * rotateY(rotation_vec.y) * rotateZ(rotation_vec.z) * scale(scale_vec.x, scale_vec.y, scale_vec.z);
	}

	Transform transform;
};

class Node {
public:
	Node(const std::string& n, const NodeType t, bool te = true) {
		name = n;
		type = t;
		transformEnabled = te;
	}

	virtual void updateBounds() {}
	virtual void setTransform(Transform t) {
		transform_data.transform = t;
		updateBounds();
	}

	virtual Transform getWorldTransform() {
		if (parent != nullptr) {
			return parent->getWorldTransform() * transform_data.transform;
		}

		return transform_data.transform;
	}

	virtual Bounds3f getWorldBounds() {
		return transform_data.transform(bounds);
	}

	bool isTransformEnabled() { return transformEnabled; }

	std::string name;
	NodeType type;

	Bounds3f bounds;
	TransformData transform_data;
	Node* parent = nullptr;

private:
	bool transformEnabled;
};

} // namespace civet

#endif // CIVET_NODE_H
