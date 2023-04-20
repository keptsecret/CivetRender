#ifndef CIVET_SCENE_H
#define CIVET_SCENE_H

#include <core/civet.h>
#include <core/mesh.h>
#include <core/light.h>

namespace civet {

class Scene {
public:
	Scene() {}

	Scene(const char* path);

	void loadScene(const char* path);
	void drawSceneTree();

	void addNode(std::shared_ptr<Node> node, NodeType type);

public:
	std::vector<std::shared_ptr<GLModel>> models;

	std::vector<std::shared_ptr<GLDirectionalLight>> dir_lights;
	std::vector<std::shared_ptr<GLPointLight>> point_lights;

	std::vector<std::shared_ptr<Node>> nodes;

private:
	void drawTreeChildren(ImGuiTreeNodeFlags node_flags, bool node_open, int index);

	bool showSceneTree = true;
};

} // namespace civet

#endif // CIVET_SCENE_H
