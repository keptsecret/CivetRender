#include <cameras/perspective.h>

namespace civet {

PerspectiveCamera::PerspectiveCamera(const AnimatedTransform& ctw, const Bounds2f& screen_window,
		float open, float close, float lens_radius, float focal_distance, float fov, Film* f, const Medium* m) :
		ProjectiveCamera(ctw, perspective(fov, 1e-2f, 1000.0f), screen_window, open, close, lens_radius, focal_distance, f, m) {
	dx_camera = raster_to_camera(Point3f(1, 0, 0)) - raster_to_camera(Point3f(0, 0, 0));
	dy_camera = raster_to_camera(Point3f(0, 1, 0)) - raster_to_camera(Point3f(0, 0, 0));
}

float PerspectiveCamera::generateRay(const CameraSample& sample, Ray* ray) const {
	Point3f p_film = Point3f(sample.p_film.x, sample.p_film.y, 0);
	Point3f p_camera = raster_to_camera(p_film);

	*ray = Ray(Point3f(0, 0, 0), normalize(Vector3f(p_camera)));
	///< modify ray for depth of field
	if (lens_radius > 0) {
		Point2f p_lens = lens_radius * concentricSamplingDisk(sample.p_lens);

		float ft = focal_distance / ray->d.z;
		Point3f p_focus = (*ray)(ft);

		ray->o = Point3f(p_lens.x, p_lens.y, 0);
		ray->d = normalize(p_focus - ray->o);
	}

	ray->time = lerp(sample.time, shutter_open, shutter_close);
	ray->medium = medium;
	*ray = camera_to_world(*ray);
	return 1;
}
float PerspectiveCamera::generateRayDifferential(const CameraSample& sample, RayDifferential* rd) const {
	Point3f p_film = Point3f(sample.p_film.x, sample.p_film.y, 0);
	Point3f p_camera = raster_to_camera(p_film);

	*rd = RayDifferential(Point3f(0, 0, 0), normalize(Vector3f(p_camera)));
	if (lens_radius > 0) {
		Point2f p_lens = lens_radius * concentricSamplingDisk(sample.p_lens);

		float ft = focal_distance / rd->d.z;
		Point3f p_focus = (*rd)(ft);

		rd->o = Point3f(p_lens.x, p_lens.y, 0);
		rd->d = normalize(p_focus - rd->o);
	}

	if (lens_radius > 0) {
		// account for lens shape and dof in differentials
		Point2f p_lens = lens_radius * concentricSamplingDisk(sample.p_lens);

		float ft = focal_distance / rd->d.z;

		Point3f p_focus = p_camera + dx_camera + (ft * Vector3f(0, 0, 1));
		rd->rx_origin = Point3f(p_lens.x, p_lens.y, 0);
		rd->rx_direction = normalize(p_focus - rd->rx_origin);

		p_focus = p_camera + dy_camera + (ft * Vector3f(0, 0, 1));
		rd->ry_origin = Point3f(p_lens.x, p_lens.y, 0);
		rd->ry_direction = normalize(p_focus - rd->ry_origin);
	} else {
		rd->rx_origin = rd->ry_origin = rd->o;
		rd->rx_direction = normalize(Vector3f(p_camera) + dx_camera);
		rd->ry_direction = normalize(Vector3f(p_camera) + dy_camera);
	}

	rd->time = lerp(sample.time, shutter_open, shutter_close);
	rd->has_differentials = true;
	rd->medium = medium;
	*rd = camera_to_world(*rd);
	return 1;
}

} // namespace civet