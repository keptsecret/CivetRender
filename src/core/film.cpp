#include <core/film.h>

namespace civet {

Film::Film(const Point2i& resolution, const Bounds2f& crop_window, std::unique_ptr<Filter> fi, float diag, const std::string& fn, float s) :
		full_resolution(resolution), diagonal(diag * 0.001f), filter(std::move(fi)), filename(fn), scale(s) {
	cropped_pixel_bounds = Bounds2i(Point2i(std::ceil(full_resolution.x * crop_window.p_min.x), std::ceil(full_resolution.y * crop_window.p_min.y)),
			Point2i(std::ceil(full_resolution.x * crop_window.p_max.x), std::ceil(full_resolution.y * crop_window.p_max.y)));
	pixels = std::unique_ptr<Pixel[]>(new Pixel[cropped_pixel_bounds.area()]);

	int offset = 0;
	for (int y = 0; y < filter_table_width; y++) {
		for (int x = 0; x < filter_table_width; x++, offset++) {
			Point2f p;
			p.x = (x + 0.5f) * filter->radius.x / filter_table_width;
			p.y = (y + 0.5f) * filter->radius.y / filter_table_width;
			filter_table[offset] = filter->evaluate(p);
		}
	}
}

Bounds2i Film::getSampleBounds() const {
	Bounds2f bounds(floor(Point2f(cropped_pixel_bounds.p_min) + Vector2f(0.5f, 0.5f) - filter->radius),
			ceil(Point2f(cropped_pixel_bounds.p_max) - Vector2f(0.5f, 0.5f) - filter->radius));
	return Bounds2i(bounds);
}

Bounds2f Film::getPhysicalExtent() const {
	///< mainly used by unimplemented realistic camera
	float aspect = float(full_resolution.y) / float(full_resolution.x);
	float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
	float y = aspect * x;
	return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

std::unique_ptr<FilmTile> Film::getFilmTile(const Bounds2i& sample_bounds) {
	Vector2f half_pixel = Vector2f(0.5f, 0.5f);
	Bounds2f bounds = Bounds2f(sample_bounds);
	Point2i p0 = (Point2i)ceil(bounds.p_min - half_pixel - filter->radius);
	Point2i p1 = (Point2i)floor(bounds.p_max - half_pixel + filter->radius) + Point2i(1, 1);
	Bounds2i tile_pixel_bounds = bIntersect(Bounds2i(p0, p1), cropped_pixel_bounds);
	return std::unique_ptr<FilmTile>(new FilmTile(tile_pixel_bounds, filter->radius, filter_table, filter_table_width));
}

void Film::mergeFilmTile(std::unique_ptr<FilmTile> tile) {
	std::lock_guard<std::mutex> lock(mutex);
	for (auto pixel : tile->getPixelBounds()) {
		// merge
		const FilmTilePixel& tile_pixel = tile->getPixel(pixel);
		Pixel& merge_pixel = getPixel(pixel);
		float xyz[3];
		tile_pixel.contrib_sum.toXYZ(xyz);
		for (int i = 0; i < 3; i++) {
			merge_pixel.xyz[i] += xyz[i];
		}
		merge_pixel.filter_weight_sum += tile_pixel.filter_weight_sum;
	}
}

void Film::setImage(const Spectrum* img) const {
	int n_pixels = cropped_pixel_bounds.area();
	for (int i = 0; i < n_pixels; i++) {
		Pixel& p = pixels[i];
		img[i].toXYZ(p.xyz);
		p.filter_weight_sum = 1;
		p.splat_xyz[0] = p.splat_xyz[1] = p.splat_xyz[2] = 0;
	}
}

void Film::addSplat(const Point2f& p, const Spectrum& v) {
	if (!bInsideExclusive(Point2i(p), cropped_pixel_bounds)) {
		return;
	}
	float xyz[3];
	v.toXYZ(xyz);
	Pixel& pixel = getPixel(Point2i(p));
	for (int i = 0; i < 3; i++) {
		pixel.splat_xyz[i].add(xyz[i]);
	}
}

void Film::writeImage(float splat_scale) {
	std::unique_ptr<float[]> rgb(new float[3 * cropped_pixel_bounds.area()]);
	int offset = 0;
	for (auto p : cropped_pixel_bounds) {
		// convert pixel XYZ to RGB
		Pixel& pixel = getPixel(p);
		XYZtoRGB(pixel.xyz, &rgb[3 * offset]);

		// normalize pixel weight sum
		float filter_wt_sum = pixel.filter_weight_sum;
		if (filter_wt_sum != 0) {
			float inv_wt = 1.0f / filter_wt_sum;
			rgb[3 * offset] = std::max(0.0f, rgb[3 * offset] * inv_wt);
			rgb[3 * offset + 1] = std::max(0.0f, rgb[3 * offset + 1] * inv_wt);
			rgb[3 * offset + 2] = std::max(0.0f, rgb[3 * offset + 2] * inv_wt);
		}

		// add splat value
		float splat_rgb[3];
		float splat_xyz[3] = { pixel.splat_xyz[0], pixel.splat_xyz[1], pixel.splat_xyz[2] };
		XYZtoRGB(splat_xyz, splat_rgb);
		rgb[3 * offset] += splat_scale * splat_rgb[0];
		rgb[3 * offset + 1] += splat_scale * splat_rgb[1];
		rgb[3 * offset + 2] += splat_scale * splat_rgb[2];

		// scale pixel value
		rgb[3 * offset] *= scale;
		rgb[3 * offset + 1] *= scale;
		rgb[3 * offset + 2] *= scale;

		offset++;
	}

	// TODO: call writeImage in ResourceManager when implemented
}

FilmTile::FilmTile(const Bounds2i& pb, const Vector2f& fr, const float* ft, int table_size) :
		pixel_bounds(pb), filter_radius(fr), inv_filter_radius(1 / filter_radius.x, 1 / filter_radius.y), filter_table(ft), filter_table_size(table_size) {
	pixels = std::vector<FilmTilePixel>(std::max(0, pixel_bounds.area()));
}

void FilmTile::addSample(const Point2f& p_film, const Spectrum& L, float sample_wt) {
	// compute bounds
	Point2f p_film_discrete = p_film - Vector2f(0.5f, 0.5f);
	Point2i p0 = (Point2i)ceil(p_film_discrete - filter_radius);
	Point2i p1 = (Point2i)floor(p_film_discrete + filter_radius) + Point2i(1, 1);
	p0 = max(p0, pixel_bounds.p_min);
	p1 = min(p1, pixel_bounds.p_max);

	// precompute offsets
	int* ifx = ALLOCA(int, p1.x - p0.x);
	for (int x = p0.x; x < p1.x; x++) {
		float fx = std::abs((x - p_film_discrete.x) * inv_filter_radius.x * filter_table_size);
		ifx[x - p0.x] = std::min((int)std::floor(fx), filter_table_size - 1);
	}
	int* ify = ALLOCA(int, p1.y - p0.y);
	for (int y = p0.y; y < p1.y; y++) {
		float fy = std::abs((y - p_film_discrete.y) * inv_filter_radius.y * filter_table_size);
		ifx[y - p0.y] = std::min((int)std::floor(fy), filter_table_size - 1);
	}

	for (int y = p0.y; y < p1.y; y++) {
		for (int x = p0.x; x < p1.x; x++) {
			// evaluate filter value at x,y pixel
			int offset = ify[y - p0.y] * filter_table_size + ifx[x - p0.x];
			float filter_wt = filter_table[offset];

			FilmTilePixel& pixel = getPixel(Point2i(x, y));
			pixel.contrib_sum += L * sample_wt * filter_wt;
			pixel.filter_weight_sum += filter_wt;
		}
	}
}

FilmTilePixel& FilmTile::getPixel(const Point2i& p) {
	int width = pixel_bounds.p_max.x - pixel_bounds.p_min.x;
	int offset = (p.x - pixel_bounds.p_min.x) + (p.y - pixel_bounds.p_min.y) * width;
	return pixels[offset];
}

} // namespace civet