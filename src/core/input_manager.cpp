#include <core/input_manager.h>

namespace civet {

std::vector<InputManager*> InputManager::instances;

void InputManager::init(GLFWwindow* window) {
	glfwSetKeyCallback(window, InputManager::keyCallback);
	glfwSetCursorPosCallback(window, InputManager::mouseCallback);
}

void InputManager::update() {
	for (auto im : instances) {
		im->setMouseOffset();
	}
}

bool InputManager::isKeyDown(unsigned int key_id) {
	bool result = false;
	auto it = key_map.find(key_id);
	if (it != key_map.end()) {
		result = key_map[key_id];
	}
	return result;
}

void InputManager::setKeyDown(unsigned int key_id, bool is_down) {
	key_map[key_id] = is_down;
}

void InputManager::setMouseCoords(Vector2f coords) {
	if (first_mouse_input) {
		mouse_coords = coords;
		first_mouse_input = false;
	}

	mouse_coords = coords;
}

void InputManager::setMouseOffset() {
	mouse_offset = mouse_coords - prev_mouse_coords;
	prev_mouse_coords = mouse_coords;

	if (!invert_y) {
		mouse_offset.y = -mouse_offset.y;
	}
}

void InputManager::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	for (auto im : instances) {
		im->setKeyDown(key, action != GLFW_RELEASE);
	}
}

void InputManager::mouseCallback(GLFWwindow* window, double x_pos_in, double y_pos_in) {
	for (auto im : instances) {
		im->setMouseCoords(Vector2f(x_pos_in, y_pos_in));
	}
}

} // namespace civet