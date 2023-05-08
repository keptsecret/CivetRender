#ifndef CIVET_EDITOR_H
#define CIVET_EDITOR_H

#include <core/civet.h>
#include <core/scene.h>

namespace civet {

class Editor {
protected:
	Editor() {}

public:
	static Editor* getSingleton();

	void draw(Scene& active_scene);

private:
	void sceneTree(Scene& active_scene);
	void sceneTreeNode(Scene& active_scene, ImGuiTreeNodeFlags node_flags, bool node_open, int index);
	void inspector(Scene& active_scene);

	bool show_scene_tree = true;
	bool show_inspector = false;
};

} // namespace civet

#endif // CIVET_EDITOR_H
