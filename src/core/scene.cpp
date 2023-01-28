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

	models.emplace_back(scene, path);

	for (int i = 0; i < scene->mNumLights; i++) {
		aiLight* light = scene->mLights[i];

		switch (light->mType) {
			case aiLightSource_DIRECTIONAL: {
				GLDirectionalLight new_light(Vector3f(light->mDirection.x, light->mDirection.y, light->mDirection.z));
				new_light.color = Vector3f(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);
				new_light.init();
				dir_lights.push_back(new_light);
				break;
			}
			case aiLightSource_POINT: {
				GLPointLight new_light(Point3f (light->mPosition.x, light->mPosition.y, light->mPosition.z));
				new_light.color = Vector3f(light->mColorDiffuse.r, light->mColorDiffuse.g, light->mColorDiffuse.b);
				new_light.attenuation.constant = light->mAttenuationConstant;
				new_light.attenuation.linear = light->mAttenuationLinear;
				new_light.attenuation.quadratic = light->mAttenuationQuadratic;
				new_light.init();
				point_lights.push_back(new_light);
				break;
			}
			default:
				break;
		}
	}
}

} // namespace civet