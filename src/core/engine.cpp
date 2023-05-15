#include <core/engine.h>

#include <stb/stb_image.h>

namespace civet {

void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window, GLCamera& camera, float delta_time) {
	InputManager::update();

	auto& io = ImGui::GetIO();
	if (io.WantCaptureMouse || io.WantCaptureKeyboard) {
		return;
	}

	Engine* engine = Engine::getSingleton();
	if (engine->input_manager.isKeyPressed(GLFW_KEY_ESCAPE)) {
		engine->active_scene.clearSelectedNode();
	}

	if (engine->input_manager.isButtonDown(GLFW_MOUSE_BUTTON_2)) { // hold right click to interact

		if (engine->input_manager.isKeyDown(GLFW_KEY_W)) {
			camera.translate(FORWARD, delta_time);
		}
		if (engine->input_manager.isKeyDown(GLFW_KEY_S)) {
			camera.translate(BACKWARD, delta_time);
		}
		if (engine->input_manager.isKeyDown(GLFW_KEY_A)) {
			camera.translate(LEFT, delta_time);
		}
		if (engine->input_manager.isKeyDown(GLFW_KEY_D)) {
			camera.translate(RIGHT, delta_time);
		}

		Vector2f mouse_offset = engine->input_manager.getMouseOffset();
		engine->view_camera.pan(mouse_offset.x, mouse_offset.y);
	}

	if (engine->input_manager.isKeyPressed(GLFW_KEY_H)) {
		Editor* editor = Editor::getSingleton();
		editor->toggleShowEditor();
	}
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

	// glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	stbi_set_flip_vertically_on_load(false);

	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_MULTISAMPLE);
	glCheckError("ERROR::Engine::init: OpenGL error code");

	DeferredRenderer::getSingleton()->init(width, height);
	parallelInit();

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	(void)io;

	ImGui::StyleColorsDark();
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init("#version 420");

	return 0;
}

int Engine::start() {
	/**
	 * Test scene
	 */
	active_scene = Scene("../civet/resources/sponza/sponza.obj");
	active_scene.models[0]->transform_data.scale_vec = Vector3f{ 0.1, 0.1, 0.1 };
	active_scene.models[0]->transform_data.updateTransform();

	const unsigned int SHADOW_RES = 4096;
	auto dir_light = std::make_shared<GLDirectionalLight>("Light_1", Vector3f(2, -1, -2), SHADOW_RES);
	dir_light->init();
	active_scene.addNode(dir_light, DirectionalLight);

	auto point_light = std::make_shared<GLPointLight>("Light_2", Point3f(0, 1, 2), SHADOW_RES);
	point_light->init();
	active_scene.addNode(point_light, PointLight);

	DeferredRenderer* renderer = DeferredRenderer::getSingleton();
	Editor* editor = Editor::getSingleton();

	renderer->setCamera(&view_camera);

	float last_frame = 0.0f;

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		renderer->setViewMat(view_camera.getViewTransform());
		Transform projection = perspective(view_camera.zoom, width / height, view_camera.near_plane, view_camera.far_plane);
		renderer->setProjectionMat(projection);

		renderer->draw(active_scene);
		editor->draw(active_scene);

		glCheckError("After draw");

		glfwPollEvents();
		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	parallelCleanup();

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

Engine* Engine::getSingleton() {
	static Engine singleton;
	return &singleton;
}

} // namespace civet