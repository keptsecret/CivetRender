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

	if (ImGui::TreeNode("Scene")) {
		ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 3);
		for (int i = 0; i < active_scene.nodes.size(); i++) {
			ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
			// only models have children for now
			if (active_scene.nodes[i]->type != Model) {
				node_flags = ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
				ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "%s", active_scene.nodes[i]->name.c_str());
			} else {
				bool node_open = ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "%s", active_scene.nodes[i]->name.c_str());

				sceneTreeNode(active_scene, node_flags, node_open, i);
			}
			if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup) && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
				active_scene.selected_node = active_scene.nodes[i];
				show_inspector = true;
			}
		}
		ImGui::TreePop();
		ImGui::PopStyleVar();
	}

	ImGui::End();
}

void Editor::sceneTreeNode(Scene& active_scene, ImGuiTreeNodeFlags node_flags, bool node_open, int index) {
	if (node_open) {
		auto n = active_scene.nodes[index];
		switch (n->type) {
			case Model:
				auto m = std::static_pointer_cast<GLModel>(n);
				auto meshes = m->getMeshes();
				for (int i = 0; i < meshes.size(); i++) {
					// meshes are leaf nodes, for now
					node_flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
					ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "%s", meshes[i].name.c_str());
				}
				break;
		}
		ImGui::TreePop();
	}
}

void Editor::inspector(Scene& active_scene) {
	if (show_inspector) {
		ImGui::Begin("Inspector", &show_inspector, ImGuiWindowFlags_HorizontalScrollbar);

		ImGui::Text("Selected item: %s", active_scene.selected_node->name.c_str());

		ImGui::End();
	}
}

} // namespace civet