#ifndef CIVET_INPUT_MANAGER_H
#define CIVET_INPUT_MANAGER_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>
#include <unordered_map>
#include <vector>

namespace civet {

class InputManager {
public:
	InputManager(float window_width, float window_height, bool inverted_y = false) :
			prev_mouse_coords(window_width / 2, window_height / 2), invert_y(inverted_y) {
		InputManager::instances.push_back(this);
	}
	~InputManager() {
		instances.erase(std::remove(instances.begin(), instances.end(), this), instances.end());
	}

	static void init(GLFWwindow* window);
	static void update();

	bool isKeyDown(unsigned int key_id);

	Vector2f getMouseOffset() { return mouse_offset; }
	Vector2f getMousePosition() { return mouse_coords; }

private:
	void setKeyDown(unsigned int key_id, bool is_down);
	//	bool wasKeyDown(unsigned int key_id);
	void setMouseCoords(Vector2f coords);
	void setMouseOffset();

	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void mouseCallback(GLFWwindow* window, double x_pos_in, double y_pos_in);

	std::unordered_map<unsigned int, bool> key_map;
	std::unordered_map<unsigned int, bool> prev_key_map;

	Vector2f mouse_coords;
	Vector2f prev_mouse_coords;
	Vector2f mouse_offset;
	bool invert_y;
	bool first_mouse_input = true;

	static std::vector<InputManager*> instances;
};

} // namespace civet

#endif // CIVET_INPUT_MANAGER_H
