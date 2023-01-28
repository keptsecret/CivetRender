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
	active_scene = Scene("../civet/resources/backpack/backpack.obj");

	const unsigned int SHADOW_RES = 4096;
	active_scene.dir_lights.push_back(GLDirectionalLight(Vector3f(2, -1, -2), SHADOW_RES));
	active_scene.dir_lights[0].init();

	active_scene.point_lights.push_back(GLPointLight(Point3f(0, 1, 2), SHADOW_RES));
	active_scene.point_lights[0].init();

	DeferredRenderer* renderer = DeferredRenderer::getSingleton();

	renderer->setCamera(&view_camera);

	float last_frame = 0.0f;

	bool show_demo_window = true;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		// TODO: implement as editor class + subclasses
		{
			ImGui_ImplOpenGL3_NewFrame();
			ImGui_ImplGlfw_NewFrame();
			ImGui::NewFrame();

			if (show_demo_window)
				ImGui::ShowDemoWindow(&show_demo_window);

			{
				static float f = 0.0f;
				static int counter = 0;

				ImGui::Begin("Hello, world!"); // Create a window called "Hello, world!" and append into it.

				ImGui::Text("This is some useful text."); // Display some text (you can use a format strings too)
				ImGui::Checkbox("Demo Window", &show_demo_window); // Edit bools storing our window open/close state

				ImGui::SliderFloat("float", &f, 0.0f, 1.0f); // Edit 1 float using a slider from 0.0f to 1.0f
				ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

				if (ImGui::Button("Button")) // Buttons return true when clicked (most widgets return true when edited/activated)
					counter++;
				ImGui::SameLine();
				ImGui::Text("counter = %d", counter);

				ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
				ImGui::End();
			}
			ImGui::Render();
		}

		Transform model = scale(0.5f, 0.5f, 0.5f);
		renderer->setModelMat(model);
		renderer->setViewMat(view_camera.getViewTransform());
		Transform projection = perspective(view_camera.zoom, width / height, 1e-2f, 1000.0f);
		renderer->setProjectionMat(projection);

		renderer->draw(active_scene);
		//renderer.draw(test_model, dir_lights, point_lights, &shader, SHADOW_RES, &depth_shader, &depth_cube_shader);

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData()); ///< here as well

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