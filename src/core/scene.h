#ifndef CIVET_SCENE_H
#define CIVET_SCENE_H

#include <core/civet.h>
#include <core/light.h>
#include <core/mesh.h>

#include <core/skybox.h>

namespace civet {

class Editor;

class Scene {
public:
	Scene() {}

	Scene(const char* path);

	void loadScene(const char* path);
	void addNode(std::shared_ptr<Node> node, NodeType type);

	void clearSelectedNode() { selected_node = nullptr; }

public:
	std::vector<std::shared_ptr<GLModel>> models;
	std::vector<std::shared_ptr<GLDirectionalLight>> dir_lights;
	std::vector<std::shared_ptr<GLPointLight>> point_lights;

	std::vector<std::shared_ptr<Node>> nodes;

	std::shared_ptr<Skybox> skybox;

private:
	friend Editor;

	std::shared_ptr<Node> selected_node = nullptr;
};

} // namespace civet

#endif // CIVET_SCENE_H
