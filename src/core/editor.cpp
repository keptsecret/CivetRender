#include <core/editor.h>

namespace civet {

Editor* Editor::getSingleton() {
	static Editor singleton;
	return &singleton;
}

void Editor::draw(Scene& active_scene) {
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	ImGui::SetNextWindowSize(ImVec2(200, 500));
	sceneTree(active_scene);

	ImGui::SetNextWindowSize(ImVec2(300, 400));
	inspector(active_scene);

	//	ImGui::Begin("Debug");
	//	ImGui::Text("Camera pitch: %.3f, yaw: %.3f", view_camera.pitch, view_camera.yaw);
	//	ImGui::Text("Camera front x: %.3f, y: %.3f, z: %.3f", view_camera.front.x, view_camera.front.y, view_camera.front.z);
	//	ImGui::Text("Camera up x: %.3f, y: %.3f, z: %.3f", view_camera.up.x, view_camera.up.y, view_camera.up.z);
	//	ImGui::End();

	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData()); ///< here as well
}

void Editor::sceneTree(Scene& active_scene) {
	ImGui::Begin("Scene Tree", &show_scene_tree, ImGuiWindowFlags_HorizontalScrollbar);

	if (ImGui::TreeNodeEx("Scene", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 1.5);
		ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2{ 0.0f, 0.0f });
		ImGui::PushStyleVar(ImGuiStyleVar_ItemInnerSpacing, ImVec2{ 3.0f, 3.0f });

		for (int i = 0; i < active_scene.nodes.size(); i++) {
			auto node = active_scene.nodes[i];
			TreeNodeState state = sceneTreeNode(active_scene, node);

			if (state.is_open) {
				if (node->type == Model) {
					auto model = std::static_pointer_cast<GLModel>(node);
					auto meshes = model->getMeshes();
					for (int j = 0; j < meshes.size(); j++) {
						sceneTreeNode(active_scene, meshes[j]);
					}
				}

				if (state.should_pop) {
					ImGui::TreePop();
				}
			}
		}

		ImGui::TreePop();
		ImGui::PopStyleVar(3);
	}

	ImGui::End();
}

TreeNodeState Editor::sceneTreeNode(Scene& active_scene, std::shared_ptr<Node> node) {
	bool is_leaf = node->type != Model;

	ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_SpanAvailWidth |
			(is_leaf
							? (ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen)
							: ImGuiTreeNodeFlags_OpenOnArrow) |
			(active_scene.selected_node == node
							? ImGuiTreeNodeFlags_Selected
							: ImGuiTreeNodeFlags_None);

	bool node_open = ImGui::TreeNodeEx(node->name.c_str(), node_flags);

	if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup) && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
		active_scene.selected_node = node;
		show_inspector = true;
	}

	return TreeNodeState{ node_open, !is_leaf };
}

void Editor::inspector(Scene& active_scene) {
	if (show_inspector) {
		ImGui::Begin("Inspector", &show_inspector, ImGuiWindowFlags_HorizontalScrollbar);

		ImGui::Text("Selected item: %s", active_scene.selected_node->name.c_str());

		ImGui::End();
	}
}

} // namespace civet