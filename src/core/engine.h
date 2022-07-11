#ifndef CIVET_ENGINE_H
#define CIVET_ENGINE_H

#include <core/civet.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

namespace civet {

class Engine {
protected:
	Engine() {}

public:
	static Engine& getSingleton();

	int init();
	int start();

	GLFWwindow* window = nullptr;
};

} // namespace civet

#endif // CIVET_ENGINE_H
