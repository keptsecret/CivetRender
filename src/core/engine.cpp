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
	Shader depth_shader("../civet/src/shaders/light_depth_vert.glsl", "../civet/src/shaders/light_depth_frag.glsl");
	Shader depth_cube_shader("../civet/src/shaders/light_cube_depth_vert.glsl", "../civet/src/shaders/light_cube_depth_frag.glsl", "../civet/src/shaders/light_cube_depth_geom.glsl");
	GLModel test_model("../civet/resources/backpack/backpack.obj");

	const unsigned int SHADOW_RES = 2048;
	std::vector<GLDirectionalLight> dir_lights;
	dir_lights.push_back(GLDirectionalLight(Vector3f(2, -1, -2), SHADOW_RES));
	dir_lights[0].init();

	std::vector<GLPointLight> point_lights;
	point_lights.push_back(GLPointLight(Point3f(0, 0, 2), SHADOW_RES));
	point_lights[0].init();

	shader.use();
	shader.setDirectionalLights(dir_lights);
	shader.setPointLights(point_lights);

	float last_frame = 0.0f;

	const unsigned int SHADOW_WIDTH = 2048, SHADOW_HEIGHT = 2048;

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Render depth map
		float near_plane = -5.0f, far_plane = 5.0f;
		Transform model = translate(Vector3f(0, 0, 0));

		depth_shader.use();
		// loop through lights and generate shadow maps
		glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
		glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
		for (auto& light : dir_lights) {
			if (light.should_update) {
				depth_shader.setMat4("model", model.m);
				light.generateShadowMap(depth_shader, near_plane, far_plane);
				test_model.draw(depth_shader, 2); ///< change tex_offset possibly
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
			}
		}

		// Render depth cube map
		float near = 1.0f, far = 25.0f;

		depth_cube_shader.use();
		// loop through lights and generate shadow maps for point lights
		glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
		glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
		for (auto& light : point_lights) {
			if (light.should_update) {
				depth_cube_shader.setMat4("model", model.m);
				light.generateShadowMap(depth_cube_shader, near, far);
				test_model.draw(depth_cube_shader, 2);
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
			}
		}

		// Render the scene
		glViewport(0, 0, width, height);
		glCullFace(GL_BACK);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use();
		Transform projection = perspective(view_camera.zoom, width / height, 1e-2f, 1000.0f);
		shader.setMat4("projection", projection.m);

		Transform view = view_camera.getViewTransform();
		shader.setMat4("view", view.m);

		shader.setMat4("model", model.m);

		shader.setVec3("viewPos", Vector3f(view_camera.position));
		shader.setFloat("material.shininess", 64.0f);

		int tex_offset = 0;
		///< put this into bind lights function in shader, possibly
		{
			for (int i = 0; i < dir_lights.size(); i++) {
				shader.setMat4(("dirLights[" + std::to_string(i) + "].light_space_mat").c_str(), dir_lights[i].light_space_mat.m);
				dir_lights[i].bindShadowMap(shader, ("dirLights[" + std::to_string(i) + "].shadow_map").c_str(), tex_offset++);
			}

			for (int i = 0; i < point_lights.size(); i++) {
				shader.setFloat(("pointLights[" + std::to_string(i) + "].far_plane").c_str(), point_lights[i].far_plane);
				point_lights[i].bindShadowMap(shader, ("pointLights[" + std::to_string(i) + "].shadow_map").c_str(), tex_offset++);
			}
		}

		test_model.draw(shader, tex_offset);
		glCheckError("After draw");

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