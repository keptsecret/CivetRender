#include <core/engine.h>

#include <stb/stb_image.h>
#include <shaders/solid.h>
#include <shaders/simple_forward.h>

namespace civet {

void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window, GLCamera& camera, float delta_time) {
	InputManager::update();

	Engine* e = Engine::getSingleton();
	if (e->input_manager.isKeyDown(GLFW_KEY_ESCAPE)) {
		glfwSetWindowShouldClose(window, true);
	}

	if (e->input_manager.isKeyDown(GLFW_KEY_W)) {
		camera.translate(FORWARD, delta_time);
	}
	if (e->input_manager.isKeyDown(GLFW_KEY_S)) {
		camera.translate(BACKWARD, delta_time);
	}
	if (e->input_manager.isKeyDown(GLFW_KEY_A)) {
		camera.translate(LEFT, delta_time);
	}
	if (e->input_manager.isKeyDown(GLFW_KEY_D)) {
		camera.translate(RIGHT, delta_time);
	}

	Vector2f mouse_offset = e->input_manager.getMouseOffset();
	e->view_camera.pan(mouse_offset.x, mouse_offset.y);
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

	InputManager::init(window);
	view_camera = GLCamera(Point3f(0, 0, 3));

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	stbi_set_flip_vertically_on_load(true);

	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_MULTISAMPLE);

	return 0;
}

int Engine::start() {
	/**
	 * Test scene
	 */
	SimpleForwardShader shader;
	GLModel test_model("../civet/resources/backpack/backpack.obj");

	std::vector<Point3f> light_pos;
	light_pos.push_back(Point3f(0, 5, 2));

	float last_frame = 0.0f;

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use(light_pos);
		Transform projection = perspective(view_camera.zoom, width / height, 1e-2f, 1000.0f);
		shader.setMat4("projection", projection.m);

		Transform view = view_camera.getViewTransform();
		shader.setMat4("view", view.m);

		Transform model = translate(Vector3f(0, 0, 0));
		shader.setMat4("model", model.m);

		shader.setVec3("viewPos", Vector3f(view_camera.position));
		shader.setFloat("material.shininess", 32.0f);
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