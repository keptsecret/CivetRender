#ifndef CIVET_ENGINE_H
#define CIVET_ENGINE_H

#include <core/camera.h>
#include <core/civet.h>
#include <core/input_manager.h>
#include <core/mesh.h>
#include <core/shader.h>

namespace civet {

class Engine {
protected:
	Engine() :
			input_manager(width, height) {}

public:
	static Engine* getSingleton();

	int init();
	int start();

	float width = 1920.0f, height = 1080.0f;

	GLFWwindow* window = nullptr;
	GLCamera view_camera;
	InputManager input_manager;
};

} // namespace civet

#endif // CIVET_ENGINE_H
