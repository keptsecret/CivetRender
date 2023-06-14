#include <core/editor.h>

#include <core/engine.h>
#include <core/skybox.h>
#include <ImGuiFileDialog/ImGuiFileDialog.h>

namespace civet {

Editor* Editor::getSingleton() {
	static Editor singleton;
	return &singleton;
}

void Editor::draw(Scene& active_scene) {
	if (show_editor) {
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		debugWindow(active_scene);

		ImGui::SetNextWindowSize(ImVec2(200, 500));
		sceneTree(active_scene);

		// ImGui::SetNextWindowSize(ImVec2(250, 400));
		inspector(active_scene);

		materialEditor(active_scene.selected_node);

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
	}
}

void Editor::debugWindow(Scene& active_scene) {
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	Engine* engine = Engine::getSingleton();

	ImGui::Begin("Debug");
	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
	scalarRangeButton(&engine->view_camera.mvmt_speed, 0.0f, 50.0f, 0xffffffffu, 0x00ffffffu, "Camera Speed", "##CS");
	scalarRangeButton(&engine->view_camera.zoom, 27.0f, 78.0f, 0xffffffffu, 0x00ffffffu, "Camera FOV", "##CFOV");
	ImGui::Text("Camera position %.3f %.3f %.3f", engine->view_camera.position.x, engine->view_camera.position.y, engine->view_camera.position.z);

	if (ImGui::TreeNodeEx("Clip range", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImGui::Indent(15.0f);
		scalarButton(&engine->view_camera.near_plane, 0xffffffffu, 0x00ffffffu, "Start:", "##CNP");
		scalarButton(&engine->view_camera.far_plane, 0xffffffffu, 0x00ffffffu, "End:", "##CFP");
		ImGui::Unindent(15.0f);
		ImGui::TreePop();
	}

	if (ImGui::Button("Open File Dialog")) {
		ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choose File", ".cpp,.h,.hpp", ".");
	}

	// display
	if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey"))
	{
		// action if OK
		if (ImGuiFileDialog::Instance()->IsOk())
		{
			std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
			std::string filePath = ImGuiFileDialog::Instance()->GetCurrentPath();
			// action
		}

		// close
		ImGuiFileDialog::Instance()->Close();
	}

	ImGui::End();
}

void Editor::sceneTree(Scene& active_scene) {
	ImGui::Begin("Scene Tree", &show_scene_tree, ImGuiWindowFlags_HorizontalScrollbar);

	if (ImGui::TreeNodeEx("Scene", ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_OpenOnArrow)) {
		ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 1.5);
		ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2{ 0.0f, 0.0f });
		ImGui::PushStyleVar(ImGuiStyleVar_ItemInnerSpacing, ImVec2{ 3.0f, 3.0f });

		if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByPopup) && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
			active_scene.selected_node = nullptr;
			active_scene.selected_self = true;
			show_material_editor = false;
		}

		for (int i = 0; i < active_scene.nodes.size(); i++) {
			auto node = active_scene.nodes[i];
			TreeNodeState state = sceneTreeNode(active_scene, node);

			if (state.is_open) {
				if (node->type == NodeType::Model) {
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
	bool is_leaf = node->type != NodeType::Model;

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
		active_scene.selected_self = false;
	}

	return TreeNodeState{ node_open, !is_leaf };
}

void Editor::inspector(Scene& active_scene) {
	show_inspector = (active_scene.selected_node != nullptr) || active_scene.selected_self;
	if (show_inspector) {
		ImGui::Begin("Inspector", &show_inspector, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_AlwaysAutoResize);

		if (active_scene.selected_self) {
			inspectScene(active_scene);
			ImGui::End();
			return;
		}

		ImGui::Text("Selected item: %s", active_scene.selected_node->name.c_str());

		auto node = active_scene.selected_node;

		if (node->isTransformEnabled()) {
			if (!ImGui::TreeNodeEx("Transform", ImGuiTreeNodeFlags_DefaultOpen)) {
				return;
			}

			ImGui::Indent(15.0f);

			ValueEditState state;
			if (ImGui::TreeNodeEx("Translation", ImGuiTreeNodeFlags_DefaultOpen)) {
				state.merge(scalarButton(&node->transform_data.translation.x, 0xff8888ffu, 0xff222266u, "X", "##T.X"));
				state.merge(scalarButton(&node->transform_data.translation.y, 0xff88ff88u, 0xff226622u, "Y", "##T.Y"));
				state.merge(scalarButton(&node->transform_data.translation.z, 0xffff8888u, 0xff662222u, "Z", "##T.Z"));
				ImGui::TreePop();
			}

			if (ImGui::TreeNodeEx("Rotation", ImGuiTreeNodeFlags_DefaultOpen)) {
				state.merge(angleButton(&node->transform_data.rotation_vec.x, 0xff8888ffu, 0xff222266u, "X", "##R.X"));
				state.merge(angleButton(&node->transform_data.rotation_vec.y, 0xff88ff88u, 0xff226622u, "Y", "##R.Y"));
				state.merge(angleButton(&node->transform_data.rotation_vec.z, 0xffff8888u, 0xff662222u, "Z", "##R.Z"));
				ImGui::TreePop();
			}

			if (ImGui::TreeNodeEx("Scale", ImGuiTreeNodeFlags_DefaultOpen)) {
				state.merge(scalarButton(&node->transform_data.scale_vec.x, 0xff8888ffu, 0xff222266u, "X", "##S.X"));
				state.merge(scalarButton(&node->transform_data.scale_vec.y, 0xff88ff88u, 0xff226622u, "Y", "##S.Y"));
				state.merge(scalarButton(&node->transform_data.scale_vec.z, 0xffff8888u, 0xff662222u, "Z", "##S.Z"));
				ImGui::TreePop();
			}

			if (state.value_changed || state.edit_finished) {
				// add more updates
				node->transform_data.updateTransform();
				node->updateWorldBounds();
			}

			ImGui::Unindent(15.0f);
			ImGui::TreePop();
		}

		if (node->type == NodeType::PointLight) {
			inspectPointLight(node);
		}

		if (node->type == NodeType::DirectionalLight) {
			inspectDirectionalLight(node);
		}

		if (node->type == NodeType::SkyBox) {
			inspectSkybox(node);
		}

		if (node->type == NodeType::Mesh) {
			if (ImGui::Button("Edit Material")) {
				show_material_editor = true;
			}
		}

		ImGui::End();
	}
}

void Editor::inspectPointLight(std::shared_ptr<Node> node) {
	if (!ImGui::TreeNodeEx("Point Light", ImGuiTreeNodeFlags_DefaultOpen)) {
		return;
	}

	ImGui::Indent(15.0f);

	auto light = std::static_pointer_cast<GLPointLight>(node);

	ImGui::Checkbox("Enabled", &light->active);

	ValueEditState node_state;
	if (ImGui::TreeNodeEx("Position", ImGuiTreeNodeFlags_DefaultOpen)) {
		node_state.merge(scalarButton(&light->position.x, 0xff8888ffu, 0xff222266u, "X", "##P.X"));
		node_state.merge(scalarButton(&light->position.y, 0xff88ff88u, 0xff226622u, "Y", "##P.Y"));
		node_state.merge(scalarButton(&light->position.z, 0xffff8888u, 0xff662222u, "Z", "##P.Z"));
		ImGui::TreePop();
	}

	colorEditVector3(&light->color, "Color");

	node_state.merge(scalarButton(&light->power, 0xffffffffu, 0x00ffffffu, "Power", "##P"));

	ImGui::Unindent(15.0f);
	ImGui::TreePop();
}

void Editor::inspectDirectionalLight(std::shared_ptr<Node> node) {
	if (!ImGui::TreeNodeEx("Directional Light", ImGuiTreeNodeFlags_DefaultOpen)) {
		return;
	}

	ImGui::Indent(15.0f);

	auto light = std::static_pointer_cast<GLDirectionalLight>(node);

	ImGui::Checkbox("Enabled", &light->active);

	ValueEditState node_state;
	if (ImGui::TreeNodeEx("Direction", ImGuiTreeNodeFlags_DefaultOpen)) {
		node_state.merge(scalarButton(&light->direction.x, 0xff8888ffu, 0xff222266u, "X", "##D.X"));
		node_state.merge(scalarButton(&light->direction.y, 0xff88ff88u, 0xff226622u, "Y", "##D.Y"));
		node_state.merge(scalarButton(&light->direction.z, 0xffff8888u, 0xff662222u, "Z", "##D.Z"));
		ImGui::TreePop();
	}

	colorEditVector3(&light->color, "Color");

	node_state.merge(scalarButton(&light->power, 0xffffffffu, 0x00ffffffu, "Power", "##P"));
	node_state.merge(scalarButton(&light->frustum_fitting_factor, 0xffffffffu, 0x00ffffffu, "Fitting Factor", "##F"));

	ImGui::Unindent(15.0f);
	ImGui::TreePop();
}

void Editor::inspectSkybox(std::shared_ptr<Node> node) {
	if (!ImGui::TreeNodeEx("Sky Box", ImGuiTreeNodeFlags_DefaultOpen)) {
		return;
	}

	ImGui::Indent(15.0f);

	auto sky = std::static_pointer_cast<Skybox>(node);

	ValueEditState node_state;
	if (ImGui::TreeNodeEx("Sun Direction", ImGuiTreeNodeFlags_DefaultOpen)) {
		node_state.merge(scalarButton(&sky->editing_params.sun_direction.x, 0xff8888ffu, 0xff222266u, "X", "##SD.X"));
		node_state.merge(scalarButton(&sky->editing_params.sun_direction.y, 0xff88ff88u, 0xff226622u, "Y", "##SD.Y"));
		node_state.merge(scalarButton(&sky->editing_params.sun_direction.z, 0xffff8888u, 0xff662222u, "Z", "##SD.Z"));
		ImGui::TreePop();
	}

	colorEditVector3(&sky->editing_params.ground_color, "Color");

	node_state.merge(scalarButton(&sky->editing_params.resolution, 0xffffffffu, 0x00ffffffu, "Resolution", "##Res"));

	if (ImGui::Button("Apply")) {
		sky->update(sky->editing_params);
	}
	ImGui::SameLine();
	if (ImGui::Button("Revert Changes")) {
		sky->resetEditingParameters();
	}

	ImGui::Unindent(15.0f);
	ImGui::TreePop();
}

void Editor::inspectScene(Scene& active_scene) {
	ImGui::Text("Scene built: %s", active_scene.isBuildForRT() ? "ready" : "pending build...");
	if (ImGui::Button("Build scene")) {
		active_scene.buildScene();
	}
	if (active_scene.isBuildForRT()) {
		ImGui::Text("Light probes baked: %s", active_scene.probe_grid->hasBakeData() ? "ready" : "pending bake...");
		if (ImGui::Button("Fit probe grid to scene")) {
			active_scene.probe_grid->fitGridToBounds(active_scene.worldBound());
		}
		if (ImGui::Button("Bake light probes")) {
			active_scene.probe_grid->initialize();
			active_scene.probe_grid->bake(active_scene);
		}

		if (ImGui::Button("test path tracer")) {
			active_scene.probe_grid->testPathtracer(active_scene);
		}

		ImGui::Text("Probe grid corner: %.3f %.3f %.3f", active_scene.probe_grid->corner_position.x, active_scene.probe_grid->corner_position.y, active_scene.probe_grid->corner_position.z);
		ImGui::Text("Probe grid dims: %d %d %d", active_scene.probe_grid->probe_grid_size.x, active_scene.probe_grid->probe_grid_size.y, active_scene.probe_grid->probe_grid_size.z);
	}
}

void Editor::materialEditor(std::shared_ptr<Node> node) {
	if (show_material_editor && node->type == NodeType::Mesh) {
		auto mesh = std::static_pointer_cast<GLMesh>(node);
		auto material = mesh->material;

		ImGui::Begin("Material Editor", &show_material_editor, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_AlwaysAutoResize);

		ImGui::Text("Selected item: %s", node->name.c_str());

		colorEditVector3(&material->albedo, "Base Color");

		scalarRangeButton(&material->metallic, 0.0f, 1.0f, 0xffffffffu, 0x00ffffffu, "Metallic", "##MM");
		scalarRangeButton(&material->roughness, 0.0f, 1.0f, 0xffffffffu, 0x00ffffffu, "Roughness", "##MR");
		scalarRangeButton(&material->ambient, 0.0f, 1.0f, 0xffffffffu, 0x00ffffffu, "Ambient Occlusion", "##MAO");

		if (ImGui::TreeNodeEx("Texture Maps", ImGuiTreeNodeFlags_OpenOnArrow)) {
			ImGui::Checkbox("Enable albedo texture", &material->use_albedo_map);
			ImGui::Checkbox("Enable metallic texture", &material->use_metallic_map);
			ImGui::Checkbox("Enable roughness texture", &material->use_roughness_map);
			ImGui::Checkbox("Enable AO texture", &material->use_ao_map);

			ImGui::Checkbox("Enable normal map", &material->use_normal_map);
			ImGui::Checkbox("Enable bump map", &material->use_bump_map);
			if (material->use_bump_map) {
				scalarRangeButton(&material->bump_scale, 0.0f, 5.0f, 0xffffffffu, 0x00ffffffu, "Bump Factor", "##MBS");
			}

			ImGui::TreePop();
		}

		ImGui::End();
	}
}

/****************************
 * Utility button and input fields
 */
static const float DRAG_MOUSE_THRESHOLD_FACTOR = 0.50f;

struct DragResult {
	bool text_edited = false;
	bool drag_edited = false;
};

DragResult customDragScalar(const char* const label,
		void* const p_data,
		const float v_speed = 1.0f,
		const float min = 0.0f,
		const float max = 0.0f,
		const ImGuiDataType data_type = ImGuiDataType_Float,
		const char* format = "%.3f") {
	ImGuiSliderFlags flags = 0;
	const float* const p_min = &min;
	const float* const p_max = &max;

	DragResult result;

	ImGuiWindow* window = ImGui::GetCurrentWindow();
	if (window->SkipItems) {
		return result;
	}

	ImGuiContext& context = *GImGui;
	const ImGuiStyle& style = context.Style;
	const ImGuiID id = window->GetID(label);
	const float w = ImGui::CalcItemWidth();

	const ImVec2 label_size = ImGui::CalcTextSize(label, nullptr, true);
	const ImRect frame_bb(window->DC.CursorPos, window->DC.CursorPos + ImVec2(w, label_size.y + style.FramePadding.y * 2.0f));
	const ImRect total_bb(frame_bb.Min, frame_bb.Max + ImVec2(label_size.x > 0.0f ? style.ItemInnerSpacing.x + label_size.x : 0.0f, 0.0f));

	const bool temp_input_allowed = (flags & ImGuiSliderFlags_NoInput) == 0;
	ImGui::ItemSize(total_bb, style.FramePadding.y);
	if (!ImGui::ItemAdd(total_bb, id, &frame_bb, temp_input_allowed ? ImGuiItemFlags_Inputable : 0)) {
		return result;
	}

	// Default format string when passing NULL
	if (format == nullptr) {
		format = ImGui::DataTypeGetInfo(data_type)->PrintFmt;
	}

	// Tabbing or CTRL-clicking on Drag turns it into an InputText
	const bool hovered = ImGui::ItemHoverable(frame_bb, id);
	bool temp_input_is_active = temp_input_allowed && ImGui::TempInputIsActive(id);
	if (!temp_input_is_active) {
		const bool input_requested_by_tabbing = temp_input_allowed && (context.LastItemData.StatusFlags & ImGuiItemStatusFlags_FocusedByTabbing) != 0;
		const bool clicked = (hovered && context.IO.MouseClicked[0]);
		const bool double_clicked = (hovered && context.IO.MouseClickedCount[0] == 2);
		if (input_requested_by_tabbing || clicked || double_clicked || context.NavActivateId == id || context.NavActivateInputId == id) {
			ImGui::SetActiveID(id, window);
			ImGui::SetFocusID(id, window);
			ImGui::FocusWindow(window);
			context.ActiveIdUsingNavDirMask = (1 << ImGuiDir_Left) | (1 << ImGuiDir_Right);
			if (temp_input_allowed) {
				if (input_requested_by_tabbing || (clicked && context.IO.KeyCtrl) || double_clicked || context.NavActivateInputId == id) {
					temp_input_is_active = true;
				}
			}
		}

		// Experimental: simple click (without moving) turns Drag into an InputText
		if (context.IO.ConfigDragClickToInputText && temp_input_allowed && !temp_input_is_active) {
			if (context.ActiveId == id && hovered && context.IO.MouseReleased[0] && !ImGui::IsMouseDragPastThreshold(0, context.IO.MouseDragThreshold * DRAG_MOUSE_THRESHOLD_FACTOR)) {
				context.NavActivateId = context.NavActivateInputId = id;
				context.NavActivateFlags = ImGuiActivateFlags_PreferInput;
				temp_input_is_active = true;
			}
		}
	}

	if (temp_input_is_active) {
		// Only clamp CTRL+Click input when ImGuiSliderFlags_AlwaysClamp is set
		const bool is_clamp_input = (flags & ImGuiSliderFlags_AlwaysClamp) != 0 && (p_min == nullptr || p_max == nullptr || ImGui::DataTypeCompare(data_type, p_min, p_max) < 0);
		result.text_edited = ImGui::TempInputScalar(frame_bb, id, label, data_type, p_data, format, is_clamp_input ? p_min : nullptr, is_clamp_input ? p_max : nullptr);
		return result;
	}

	// Draw frame
	const ImU32 frame_col = ImGui::GetColorU32(context.ActiveId == id ? ImGuiCol_FrameBgActive : hovered ? ImGuiCol_FrameBgHovered
																										 : ImGuiCol_FrameBg);
	ImGui::RenderNavHighlight(frame_bb, id);
	ImGui::RenderFrame(frame_bb.Min, frame_bb.Max, frame_col, true, style.FrameRounding);

	// Drag behavior
	result.drag_edited = ImGui::DragBehavior(id, data_type, p_data, v_speed, p_min, p_max, format, flags);
	if (result.drag_edited) {
		ImGui::MarkItemEdited(id);
	}

	// Display value using user-provided display format so user can add prefix/suffix/decorations to the value.
	char value_buf[64];
	const char* value_buf_end = value_buf + ImGui::DataTypeFormatString(value_buf, IM_ARRAYSIZE(value_buf), data_type, p_data, format);
	if (context.LogEnabled) {
		ImGui::LogSetNextTextDecoration("{", "}");
	}
	ImGui::RenderTextClipped(frame_bb.Min, frame_bb.Max, value_buf, value_buf_end, nullptr, ImVec2(0.5f, 0.5f));

	if (label_size.x > 0.0f) {
		ImGui::RenderText(ImVec2(frame_bb.Max.x + style.ItemInnerSpacing.x, frame_bb.Min.y + style.FramePadding.y), label);
	}

	return result;
}

ValueEditState Editor::scalarButton(float* value, uint32_t text_color, uint32_t background_color, const char* label, const char* imgui_label) const {
	ImGui::PushStyleColor(ImGuiCol_Text, text_color);
	ImGui::PushStyleColor(ImGuiCol_Button, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonHovered, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonActive, background_color);
	ImGui::SetNextItemWidth(12.0f);
	ImGui::Button(label);
	ImGui::SameLine();
	ImGui::PopStyleColor(4);
	ImGui::SetNextItemWidth(100.0f);

	constexpr float value_speed = 0.02f;
	const auto value_changed = customDragScalar(imgui_label, value, value_speed);
	const bool edit_ended = ImGui::IsItemDeactivatedAfterEdit();

	return ValueEditState{ value_changed.drag_edited || (value_changed.text_edited && edit_ended), edit_ended };
}

ValueEditState Editor::scalarButton(unsigned int* value, uint32_t text_color, uint32_t background_color, const char* label, const char* imgui_label) const {
	ImGui::PushStyleColor(ImGuiCol_Text, text_color);
	ImGui::PushStyleColor(ImGuiCol_Button, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonHovered, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonActive, background_color);
	ImGui::SetNextItemWidth(12.0f);
	ImGui::Button(label);
	ImGui::SameLine();
	ImGui::PopStyleColor(4);
	ImGui::SetNextItemWidth(100.0f);

	constexpr unsigned int value_speed = 1;
	const auto value_changed = customDragScalar(imgui_label, value, value_speed, 0, 0, ImGuiDataType_U32, nullptr);
	const bool edit_ended = ImGui::IsItemDeactivatedAfterEdit();

	return ValueEditState{ value_changed.drag_edited || (value_changed.text_edited && edit_ended), edit_ended };
}

ValueEditState Editor::scalarRangeButton(float* value, float value_min, float value_max, uint32_t text_color, uint32_t background_color, const char* label, const char* imgui_label) const {
	ImGui::PushStyleColor(ImGuiCol_Text, text_color);
	ImGui::PushStyleColor(ImGuiCol_Button, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonHovered, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonActive, background_color);
	ImGui::SetNextItemWidth(12.0f);
	ImGui::Button(label);
	ImGui::SameLine();
	ImGui::PopStyleColor(4);
	ImGui::SetNextItemWidth(100.0f);

	const float value_speed = (value_max - value_min) / 200.f;
	const auto value_changed = customDragScalar(imgui_label, value, value_speed, value_min, value_max);
	const bool edit_ended = ImGui::IsItemDeactivatedAfterEdit();

	return ValueEditState{ value_changed.drag_edited || (value_changed.text_edited && edit_ended), edit_ended };
}

ValueEditState Editor::angleButton(float* value, uint32_t text_color, uint32_t background_color, const char* label, const char* imgui_label) const {
	ImGui::PushStyleColor(ImGuiCol_Text, text_color);
	ImGui::PushStyleColor(ImGuiCol_Button, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonHovered, background_color);
	ImGui::PushStyleColor(ImGuiCol_ButtonActive, background_color);
	ImGui::SetNextItemWidth(12.0f);
	ImGui::Button(label);
	ImGui::SameLine();
	ImGui::PopStyleColor(4);
	ImGui::SetNextItemWidth(100.0f);

	float temp_value = *value;
	constexpr float value_speed = 1.0f;
	constexpr float value_min = -360.0f;
	constexpr float value_max = 360.0f;

	const auto drag_result = customDragScalar(
			imgui_label,
			&temp_value,
			value_speed,
			value_min,
			value_max,
			ImGuiDataType_Float,
			"%.f\xc2\xb0" // degree symbol in UTF-8
	);

	const bool edit_ended = ImGui::IsItemDeactivatedAfterEdit();
	const bool value_changed = drag_result.drag_edited || (drag_result.text_edited && edit_ended);
	if (value_changed) {
		*value = temp_value;
	}

	return ValueEditState{ value_changed, ImGui::IsItemDeactivatedAfterEdit() };
}

bool Editor::colorEditVector3(Vector3f* color, const char* imgui_label) {
	float col3[3] = {color->x, color->y, color->z};
	bool res = ImGui::ColorEdit3("Color", col3);
	color->x = col3[0], color->y = col3[1], color->z = col3[2];

	return res;
}

} // namespace civet