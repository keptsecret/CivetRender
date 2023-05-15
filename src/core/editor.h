#ifndef CIVET_EDITOR_H
#define CIVET_EDITOR_H

#include <core/civet.h>
#include <core/scene.h>

namespace civet {

struct TreeNodeState {
	bool is_open;
	bool should_pop;
};

struct ValueEditState {
	bool value_changed = false;
	bool edit_finished = false;

	void merge(const ValueEditState& state) {
		value_changed = value_changed || state.value_changed;
		edit_finished = edit_finished || state.edit_finished;
	}
};

class Editor {
protected:
	Editor() {}

public:
	static Editor* getSingleton();

	void draw(Scene& active_scene);

	void toggleShowEditor() { show_editor = !show_editor; }

private:
	void debugWindow(Scene& active_scene);

	void sceneTree(Scene& active_scene);
	TreeNodeState sceneTreeNode(Scene& active_scene, std::shared_ptr<Node> node);

	void inspector(Scene& active_scene);
	void inspectPointLight(std::shared_ptr<Node> node);
	void inspectDirectionalLight(std::shared_ptr<Node> node);
	void inspectSkybox(std::shared_ptr<Node> node);

	ValueEditState scalarButton(float* value, uint32_t text_color, uint32_t background_color,
			const char* label, const char* imgui_label) const;
	ValueEditState scalarButton(unsigned int* value, uint32_t text_color, uint32_t background_color,
			const char* label, const char* imgui_label) const;
	ValueEditState scalarRangeButton(float* value, float min, float max, uint32_t text_color, uint32_t background_color,
			const char* label, const char* imgui_label) const;
	ValueEditState angleButton(float* value, uint32_t text_color, uint32_t background_color,
			const char* label, const char* imgui_label) const;

	bool show_editor = true;
	bool show_scene_tree = true;
	bool show_inspector = false;
};

} // namespace civet

#endif // CIVET_EDITOR_H
