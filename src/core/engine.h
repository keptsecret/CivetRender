#ifndef CIVET_ENGINE_H
#define CIVET_ENGINE_H

#include <core/camera.h>
#include <core/civet.h>
#include <core/input_manager.h>
#include <core/shader.h>
#include <core/scene.h>
#include <rendering/deferred_renderer.h>

namespace civet {

class Engine {
protected:
	Engine() :
			input_manager(width, height),
			frame_time(1.f / MAX_FPS) {}

public:
	static Engine* getSingleton();

	int init();
	int start();

	float width = 1920.0f, height = 1080.0f;
	float MAX_FPS = 120.f;
	float frame_time;

	Scene active_scene;

	GLFWwindow* window = nullptr;
	GLCamera view_camera;
	InputManager input_manager;
};

} // namespace civet

#endif // CIVET_ENGINE_H
