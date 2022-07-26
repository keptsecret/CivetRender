#ifndef CIVET_ENGINE_H
#define CIVET_ENGINE_H

#include <core/civet.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <core/camera.h>
#include <core/shader.h>

namespace civet {

class Engine {
protected:
	Engine() {}

public:
	static Engine* getSingleton();

	int init();
	int start();

	GLFWwindow* window = nullptr;
	GLCamera view_camera;

	float width = 800.0f, height = 600.0f;

	float last_x = width / 2.0f;
	float last_y = height / 2.0f;
	bool first_mouse = true;
	bool invert_y = false;
};

} // namespace civet

#endif // CIVET_ENGINE_H
