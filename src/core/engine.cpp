#include <core/engine.h>

#include <shaders/simple_forward.h>
#include <shaders/solid.h>
#include <stb/stb_image.h>

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

	stbi_set_flip_vertically_on_load(true);

	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_MULTISAMPLE);
	glCheckError("ERROR::Engine::init: OpenGL error code");

	DeferredRenderer::getSingleton()->init(width, height);
	//renderer.init();

	return 0;
}

int Engine::start() {
	/**
	 * Test scene
	 */
	//SimpleForwardShader shader;
	//Shader shader("../civet/src/shaders/deferred_geometry_vert.glsl", "../civet/src/shaders/deferred_geometry_frag.glsl");
	Shader depth_shader("../civet/src/shaders/light_depth_vert.glsl", "../civet/src/shaders/light_depth_frag.glsl");
	Shader depth_cube_shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/light_cube_depth_frag.glsl", "../civet/src/shaders/light_cube_depth_geom.glsl");
	GLModel test_model("../civet/resources/backpack/backpack.obj");

	const unsigned int SHADOW_RES = 4096;
	std::vector<GLDirectionalLight> dir_lights;
	dir_lights.push_back(GLDirectionalLight(Vector3f(2, -1, -2), SHADOW_RES));
	dir_lights[0].init();

	std::vector<GLPointLight> point_lights;
	point_lights.push_back(GLPointLight(Point3f(0, 0, 2), SHADOW_RES));
	point_lights[0].init();

//	shader.use();
//	shader.setDirectionalLights(dir_lights);
//	shader.setPointLights(point_lights);

	DeferredRenderer* renderer = DeferredRenderer::getSingleton();

	renderer->setCamera(&view_camera);

	float last_frame = 0.0f;

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		Transform model = scale(0.5f, 0.5f, 0.5f);
		renderer->setModelMat(model);
		renderer->setViewMat(view_camera.getViewTransform());
		Transform projection = perspective(view_camera.zoom, width / height, 1e-2f, 1000.0f);
		renderer->setProjectionMat(projection);

		renderer->draw(test_model, dir_lights, point_lights);
		//renderer.draw(test_model, dir_lights, point_lights, &shader, SHADOW_RES, &depth_shader, &depth_cube_shader);

		glCheckError("After draw");

		glfwPollEvents();
		glfwSwapBuffers(window);
	}

	glfwTerminate();
	return 0;
}

Engine* Engine::getSingleton() {
	static Engine singleton;
	return &singleton;
}

} // namespace civet