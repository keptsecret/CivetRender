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

	std::vector<GLModel> models;

	std::vector<GLDirectionalLight> dir_lights;
	std::vector<GLPointLight> point_lights;
};

} // namespace civet

#endif // CIVET_SCENE_H
