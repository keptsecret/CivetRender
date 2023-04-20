#ifndef CIVET_NODE_H
#define CIVET_NODE_H

#include <core/civet.h>

namespace civet {

enum NodeType {
	Model,
	Mesh,
	DirectionalLight,
	PointLight
};

struct TransformData {
	// TODO: scene graph type transform data
};

class Node {
public:
	Node(const std::string& n, const NodeType t) {
		name = n;
		type = t;
	}

	std::string name;
	NodeType type;
};

} // namespace civet

#endif // CIVET_NODE_H
