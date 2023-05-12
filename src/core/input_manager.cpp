#include <core/input_manager.h>

namespace civet {

std::vector<InputManager*> InputManager::instances;

void InputManager::init(GLFWwindow* window) {
	glfwSetKeyCallback(window, InputManager::keyCallback);
	glfwSetMouseButtonCallback(window, InputManager::mouseButtonCallback);
	glfwSetCursorPosCallback(window, InputManager::mouseCallback);
}

void InputManager::update() {
	for (auto im : instances) {
		im->setMouseOffset();
		//im->key_pressed_map.clear();
	}
}

bool InputManager::isKeyDown(unsigned int key_id) {
	bool result = false;
	auto it = key_down_map.find(key_id);
	if (it != key_down_map.end()) {
		result = key_down_map[key_id];
	}
	return result;
}

bool InputManager::isKeyPressed(unsigned int key_id) {
	bool result = false;
	auto it = key_pressed_map.find(key_id);
	if (it != key_pressed_map.end()) {
		result = key_pressed_map[key_id];
		key_pressed_map[key_id] = false;
	}
	return result;
}

bool InputManager::isButtonDown(unsigned int butt_id) {
	bool result = false;
	auto it = mb_map.find(butt_id);
	if (it != mb_map.end()) {
		result = mb_map[butt_id];
	}
	return result;
}

void InputManager::setKeyDown(unsigned int key_id, bool is_down) {
	key_down_map[key_id] = is_down;
}

void InputManager::setKeyPressed(unsigned int key_id) {
	auto it = key_pressed_map.find(key_id);
	if (it != key_pressed_map.end()) {
		key_pressed_map[key_id] = !key_pressed_map[key_id];
	} else {
		key_pressed_map[key_id] = true;
	}
}

void InputManager::setButtonDown(unsigned int butt_id, bool is_down) {
	mb_map[butt_id] = is_down;
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
		if (action == GLFW_PRESS) {
			im->setKeyPressed(key);
		}
	}
}

void InputManager::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	for (auto im : instances) {
		im->setButtonDown(button, action != GLFW_RELEASE);
	}
}

void InputManager::mouseCallback(GLFWwindow* window, double x_pos_in, double y_pos_in) {
	for (auto im : instances) {
		im->setMouseCoords(Vector2f(x_pos_in, y_pos_in));
	}
}

} // namespace civet