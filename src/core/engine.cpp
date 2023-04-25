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

	Engine* e = Engine::getSingleton();
	if (e->input_manager.isButtonPressed(GLFW_MOUSE_BUTTON_2)) {	// hold right click to interact
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

	//glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	stbi_set_flip_vertically_on_load(true);

	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_MULTISAMPLE);
	glCheckError("ERROR::Engine::init: OpenGL error code");

	DeferredRenderer::getSingleton()->init(width, height);
	//renderer.init();

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

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
	active_scene.models[0]->setTransform(scale(0.1, 0.1, 0.1));

	const unsigned int SHADOW_RES = 4096;
	auto dir_light = std::make_shared<GLDirectionalLight>("Light_1", Vector3f(2, -1, -2), SHADOW_RES);
	dir_light->init();
	active_scene.addNode(dir_light, DirectionalLight);

//	auto point_light = std::make_shared<GLPointLight>("Light_2", Point3f(0, 1, 2), SHADOW_RES);
//	point_light->init();
//	active_scene.addNode(point_light, PointLight);

	DeferredRenderer* renderer = DeferredRenderer::getSingleton();

	renderer->setCamera(&view_camera);

	float last_frame = 0.0f;

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		renderer->setViewMat(view_camera.getViewTransform());
		Transform projection = perspective(view_camera.zoom, width / height, 1e-2f, 1000.0f);
		renderer->setProjectionMat(projection);

		renderer->draw(active_scene);

		// TODO: implement as editor class + subclasses
		{
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();

			ImGui::SetNextWindowSize(ImVec2(150, 500));
			active_scene.drawSceneTree();

			ImGui::Begin("Debug");
			ImGui::Text("Camera pitch: %.3f, yaw: %.3f", view_camera.pitch, view_camera.yaw);
			ImGui::Text("Camera front x: %.3f, y: %.3f, z: %.3f", view_camera.front.x, view_camera.front.y, view_camera.front.z);
			ImGui::Text("Camera up x: %.3f, y: %.3f, z: %.3f", view_camera.up.x, view_camera.up.y, view_camera.up.z);
			ImGui::End();

			ImGui::Render();
			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData()); ///< here as well
		}

		glCheckError("After draw");

		glfwPollEvents();
		glfwSwapBuffers(window);
	}

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

Engine* Engine::getSingleton() {
	static Engine singleton;
	return &singleton;
}

} // namespace civet