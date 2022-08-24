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

	std::vector<Point3f> light_pos;
	light_pos.push_back(Point3f(0, 0, 2));

	float last_frame = 0.0f;

	// Generate framebuffer and texture to store depth map (for directional light)
	unsigned int depth_map_FBO;
	glGenFramebuffers(1, &depth_map_FBO);

	const unsigned int SHADOW_WIDTH = 2048, SHADOW_HEIGHT = 2048;
	unsigned int depth_map;
	glGenTextures(1, &depth_map);
	glBindTexture(GL_TEXTURE_2D, depth_map);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	float border_color[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border_color);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, depth_map_FBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_map, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	shader.use(light_pos);
	shader.setInt("shadowMap", 0);

	// Framebuffer and texture for depth map (for one point light) -- using cube map
	unsigned int depth_cube_map_FBO;
	glGenFramebuffers(1, &depth_cube_map_FBO);

	unsigned int depth_cube_map;
	glGenTextures(1, &depth_cube_map);

	glBindTexture(GL_TEXTURE_CUBE_MAP, depth_cube_map);
	for (unsigned int i = 0; i < 6; i++) {
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, depth_cube_map_FBO);
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_cube_map, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	shader.setInt("shadowMap2", 1);

	while (!glfwWindowShouldClose(window)) {
		float current_frame = float(glfwGetTime());
		float delta_time = current_frame - last_frame;
		last_frame = current_frame;

		processInput(window, view_camera, delta_time);

		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Render depth map
		float near_plane = -5.0f, far_plane = 5.0f;
		Transform light_projection = orthographic(-10.0f, 10.0f, -10.0f, 10.0f, near_plane, far_plane);
		Transform light_view = lookAtRH(Point3f(2, -1, -1), Point3f(0, 0, 0), Vector3f(0, 1, 0));
		Transform light_space_tf = light_projection * light_view;

		Transform model = translate(Vector3f(0, 0, 0));

		depth_shader.use(light_pos);
		depth_shader.setMat4("model", model.m);
		depth_shader.setMat4("lightSpaceMatrix", light_space_tf.m);

		glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
		glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
		glBindFramebuffer(GL_FRAMEBUFFER, depth_map_FBO);
		glClear(GL_DEPTH_BUFFER_BIT);
		test_model.draw(depth_shader, 2);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// Render depth cube map
		float near = 1.0f, far = 25.0f;
		Transform point_projection = perspective(90.0f, 1.0f, near, far);
		std::vector<Transform> point_transforms;
		point_transforms.push_back(point_projection * lookAtRH(light_pos[0], light_pos[0] + Vector3f(1, 0, 0), Vector3f(0, -1, 0)));
		point_transforms.push_back(point_projection * lookAtRH(light_pos[0], light_pos[0] + Vector3f(-1, 0, 0), Vector3f(0, -1, 0)));
		point_transforms.push_back(point_projection * lookAtRH(light_pos[0], light_pos[0] + Vector3f(0, 1, 0), Vector3f(0, 0, 1)));
		point_transforms.push_back(point_projection * lookAtRH(light_pos[0], light_pos[0] + Vector3f(0, -1, 0), Vector3f(0, 0, -1)));
		point_transforms.push_back(point_projection * lookAtRH(light_pos[0], light_pos[0] + Vector3f(0, 0, 1), Vector3f(0, -1, 0)));
		point_transforms.push_back(point_projection * lookAtRH(light_pos[0], light_pos[0] + Vector3f(0, 0, -1), Vector3f(0, -1, 0)));

		depth_cube_shader.use(light_pos);
		for (unsigned int i = 0; i < 6; i++) {
			depth_cube_shader.setMat4("shadowMatrices[" + std::to_string(i) + "]", point_transforms[i].m);
		}
		depth_cube_shader.setMat4("model", model.m);
		depth_cube_shader.setFloat("far_plane", far);
		depth_cube_shader.setVec3("lightPos", Vector3f(light_pos[0]));

		glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
		glCullFace(GL_FRONT); ///< fix for peter panning shadow artifacts
		glBindFramebuffer(GL_FRAMEBUFFER, depth_cube_map_FBO);
		glClear(GL_DEPTH_BUFFER_BIT);
		test_model.draw(depth_cube_shader, 2);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// Render the scene
		glViewport(0, 0, width, height);
		glCullFace(GL_BACK);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		shader.use(light_pos);
		Transform projection = perspective(view_camera.zoom, width / height, 1e-2f, 1000.0f);
		shader.setMat4("projection", projection.m);

		Transform view = view_camera.getViewTransform();
		shader.setMat4("view", view.m);

		shader.setMat4("model", model.m);

		shader.setVec3("viewPos", Vector3f(view_camera.position));
		shader.setFloat("material.shininess", 64.0f);
		shader.setMat4("lightSpaceMatrix", light_space_tf.m);

		// temp shadows, need to move
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, depth_map);
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_CUBE_MAP, depth_cube_map);

		test_model.draw(shader, 2);

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