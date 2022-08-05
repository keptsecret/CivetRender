#include <core/engine.h>

namespace civet {

void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window, GLCamera& camera, float delta_time) {
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		camera.translate(FORWARD, delta_time);
	}
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
		camera.translate(BACKWARD, delta_time);
	}
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
		camera.translate(LEFT, delta_time);
	}
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
		camera.translate(RIGHT, delta_time);
	}
}

void mouseCallback(GLFWwindow* window, double x_pos_in, double y_pos_in) {
	float x_pos = float(x_pos_in);
	float y_pos = float(y_pos_in);

	Engine* e = Engine::getSingleton();
	if (e->first_mouse) {
		e->last_x = x_pos;
		e->last_y = y_pos;
		e->first_mouse = false;
	}

	float xoff = x_pos - e->last_x;
	float yoff = e->last_y - y_pos;

	if (e->invert_y) {
		yoff = -yoff;
	}

	e->last_x = x_pos;
	e->last_y = y_pos;
	e->view_camera.pan(xoff, yoff);
}

int Engine::init() {
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	window = glfwCreateWindow(width, height, "New Window", nullptr, nullptr);
	if (window == nullptr) {
		printf("ERROR::Engine: Failed to create window.\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		printf("ERROR::Engine: Failed to load GLAD.\n");
		return -1;
	}
	glViewport(0, 0, width, height);

	glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
	glfwSetCursorPosCallback(window, mouseCallback);

	view_camera = GLCamera(Point3f(0, 0, 3));

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_FRAMEBUFFER_SRGB);	///< enables automatic gamma correction

	return 0;
}

int Engine::start() {
	/**
	 * Test scene
	 */
	Shader shader("../civet/src/material/matte_shader.vert", "../civet/src/material/matte_shader.frag");

	GLModel test_model("../civet/resources/backpack/backpack.obj");

	float last_frame = 0.0f;

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use();
		Transform projection = perspective(view_camera.zoom, width / height, 0.1f, 100.0f);
		shader.setMat4("projection", projection.m);

		Transform view = view_camera.getViewTransform();
		shader.setMat4("view", view.m);

		Transform model = translate(Vector3f(0, 0, 0));
		shader.setMat4("model", model.m);

		shader.setVec3("viewPos", Vector3f(view_camera.position));
		test_model.draw(shader);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}

Engine* Engine::getSingleton() {
	static Engine singleton;
	return &singleton;
}

} // namespace civet