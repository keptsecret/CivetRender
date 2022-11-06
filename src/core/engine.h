#ifndef CIVET_ENGINE_H
#define CIVET_ENGINE_H

#include <core/camera.h>
#include <core/civet.h>
#include <core/input_manager.h>
#include <core/mesh.h>
#include <core/shader.h>
#include <rendering/forward_renderer.h>

namespace civet {

class Engine {
protected:
	Engine() :
			input_manager(width, height),
			renderer(width, height) {}

public:
	static Engine* getSingleton();

	int init();
	int start();

	float width = 1920.0f, height = 1080.0f;

	GLFWwindow* window = nullptr;
	GLCamera view_camera;
	InputManager input_manager;
	ForwardRenderer renderer;

private:
	void renderGeometry();

	void renderLights();
};

} // namespace civet

#endif // CIVET_ENGINE_H
