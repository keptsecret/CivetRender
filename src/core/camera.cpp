#include <core/camera.h>

namespace civet {

Camera::Camera(const AnimatedTransform& ctw, float open, float close, Film* f, const Medium* m) :
		camera_to_world(ctw), shutter_open(open), shutter_close(close), film(f), medium(m) {
}

float Camera::generateRayDifferential(const CameraSample& sample, RayDifferential* rd) const {
	float wt = generateRay(sample, rd);
	if (wt == 0) {
		return 0;
	}

	float wtx;
	for (float eps : { 0.05f, -0.05f }) {
		CameraSample sshift = sample;
		sshift.p_film.x += eps;
		Ray rx;
		wtx = generateRay(sshift, &rx);
		rd->rx_origin = rd->o + (rx.o - rd->o) / eps;
		rd->rx_direction = rd->d + (rx.d - rd->d) / eps;
		if (wtx != 0) {
			break;
		}
	}
	if (wtx == 0) {
		return 0;
	}

	float wty;
	for (float eps : { 0.05f, -0.05f }) {
		CameraSample sshift = sample;
		sshift.p_film.y += eps;
		Ray ry;
		wty = generateRay(sshift, &ry);
		rd->ry_origin = rd->o + (ry.o - rd->o) / eps;
		rd->ry_direction = rd->d + (ry.d - rd->d) / eps;
		if (wty != 0) {
			break;
		}
	}
	if (wty == 0) {
		return 0;
	}

	rd->has_differentials = true;
	return wt;
}

ProjectiveCamera::ProjectiveCamera(const AnimatedTransform& ctw, const Transform& cts, const Bounds2f& screen_window,
		float open, float close, float lensr, float focald, Film* f, const Medium* m) :
		Camera(ctw, open, close, f, m), camera_to_screen(cts), lens_radius(lensr), focal_distance(focald) {
	screen_to_raster = scale(film->full_resolution.x, film->full_resolution.y, 1) *
			scale(1 / (screen_window.p_max.x - screen_window.p_min.x), 1 / (screen_window.p_max.y - screen_window.p_min.y), 1) *
			translate(Vector3f(-screen_window.p_min.x, -screen_window.p_max.y, 0));
	raster_to_screen = inverse(screen_to_raster);
	raster_to_camera = inverse(camera_to_screen) * raster_to_screen;
}

GLCamera::GLCamera(Point3f pos, Vector3f up, float _yaw, float _pitch) :
		position(pos), front(Vector3f(0, 0, -1)), world_up(up), yaw(_yaw), pitch(_pitch), mvmt_speed(2.5f), sensitivity(0.1f), zoom(45.0f) {
	updateCameraVectors();
}

Transform GLCamera::getViewTransform() const {
	return lookAtRH(position, position + front, up);
}

void GLCamera::translate(CameraMovement direction, float delta_time) {
	float velocity = mvmt_speed * delta_time;
	switch (direction) {
		case FORWARD:
			position += front * velocity;
			break;
		case BACKWARD:
			position -= front * velocity;
			break;
		case LEFT:
			position -= right * velocity;
			break;
		case RIGHT:
			position += right * velocity;
			break;
	}
}

void GLCamera::pan(float xoffset, float yoffset, bool constrain_pitch) {
	xoffset *= sensitivity;
	yoffset *= sensitivity;

	yaw += xoffset;
	pitch += yoffset;

	if (constrain_pitch) {
		pitch = clamp(pitch, -89.0f, 89.0f);
	}

	updateCameraVectors();
}

void GLCamera::updateCameraVectors() {
	Vector3f direction;
	direction.x = std::cos(radians(yaw)) * std::cos(radians(pitch));
	direction.y = std::sin(radians(pitch));
	direction.z = std::sin(radians(yaw)) * std::cos(radians(pitch));
	front = normalize(direction);
	right = normalize(cross(front, world_up));
	up = normalize(cross(right, front));
}

} // namespace civet