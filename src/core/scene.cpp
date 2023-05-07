#include <core/scene.h>

#include <assimp/postprocess.h>
#include <assimp/Importer.hpp>

namespace civet {

Scene::Scene(const char* path) {
	loadScene(path);
}

void Scene::loadScene(const char* path) {
	Assimp::Importer importer;
	const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_GenNormals | aiProcess_CalcTangentSpace | aiProcess_FlipUVs);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
		std::cout << "ERROR::GLModel: assimp error: " << importer.GetErrorString() << '\n';
		return;
	}

	auto model = std::make_shared<GLModel>("Model");
	model->loadModel(scene, path);
	models.push_back(model);
	nodes.push_back(model);

	for (int i = 0; i < scene->mNumLights; i++) {
		aiLight* light = scene->mLights[i];
		std::string light_name = ("Light_" + std::to_string(i));

		switch (light->mType) {
			case aiLightSource_DIRECTIONAL: {
				auto new_light = std::make_shared<GLDirectionalLight>(light_name, Vector3f(light->mDirection.x, light->mDirection.y, light->mDirection.z));
				new_light->color = Vector3f(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);
				new_light->init();
				dir_lights.push_back(new_light);
				nodes.push_back(new_light);
				break;
			}
			case aiLightSource_POINT: {
				auto new_light = std::make_shared<GLPointLight>(light_name, Point3f(light->mPosition.x, light->mPosition.y, light->mPosition.z));
				new_light->color = Vector3f(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);
				new_light->attenuation.constant = light->mAttenuationConstant;
				new_light->attenuation.linear = light->mAttenuationLinear;
				new_light->attenuation.quadratic = light->mAttenuationQuadratic;
				new_light->init();
				point_lights.push_back(new_light);
				nodes.push_back(new_light);
				break;
			}
			default:
				break;
		}
	}
}

void Scene::addNode(std::shared_ptr<Node> node, NodeType type) {
	switch (type) {
		case NodeType::Model:
			models.push_back(std::static_pointer_cast<GLModel>(node));
			break;
		case NodeType::DirectionalLight:
			dir_lights.push_back(std::static_pointer_cast<GLDirectionalLight>(node));
			break;
		case NodeType::PointLight:
			point_lights.push_back(std::static_pointer_cast<GLPointLight>(node));
			break;
	}

	nodes.push_back(node);
}

void Scene::drawSceneTree() {
	ImGui::Begin("Scene Tree", &showSceneTree, ImGuiWindowFlags_HorizontalScrollbar);

	if (ImGui::TreeNode("Scene")) {
		ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 3);
		for (int i = 0; i < nodes.size(); i++) {
			ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
			// only models have children for now
			if (nodes[i]->type != Model) {
				node_flags = ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
				ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "%s", nodes[i]->name.c_str());
			} else {
				bool node_open = ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "%s", nodes[i]->name.c_str());

				drawTreeChildren(node_flags, node_open, i);
			}
		}
		ImGui::TreePop();
		ImGui::PopStyleVar();
	}

	ImGui::End();
}

void Scene::drawTreeChildren(ImGuiTreeNodeFlags node_flags, bool node_open, int index) {
	if (node_open) {
		auto n = nodes[index];
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

} // namespace civet