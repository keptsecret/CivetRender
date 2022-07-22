#ifndef CIVET_FILM_H
#define CIVET_FILM_H

#include <core/civet.h>
#include <core/filter.h>
#include <core/geometry/vecmath.h>
#include <core/spectrum.h>
#include <utils/parallel.h>

namespace civet {

struct FilmTilePixel {
	Spectrum contrib_sum = 0.0f;
	float filter_weight_sum = 0.0f;
};

class Film {
public:
	Film(const Point2i& resolution, const Bounds2f& crop_window, std::unique_ptr<Filter> fi, float diag, const std::string& fn, float s);

	Bounds2i getSampleBounds() const;
	Bounds2f getPhysicalExtent() const;

	std::unique_ptr<FilmTile> getFilmTile(const Bounds2i& sample_bounds);

	void mergeFilmTile(std::unique_ptr<FilmTile> tile);

	void setImage(const Spectrum* img) const;

	void addSplat(const Point2f& p, const Spectrum& v);

	void writeImage(float splat_scale);

	const Point2i full_resolution;
	const float diagonal;
	std::unique_ptr<Filter> filter;
	const std::string filename;
	Bounds2i cropped_pixel_bounds;

private:
	struct Pixel {
		float xyz[3] = { 0, 0, 0 };
		float filter_weight_sum = 0;
		AtomicFloat splat_xyz[3];
		float pad;
	};
	std::unique_ptr<Pixel[]> pixels;
	static constexpr int filter_table_width = 16;
	float filter_table[filter_table_width * filter_table_width];
	std::mutex mutex;
	const float scale;

	Pixel& getPixel(const Point2i& p) {
		int width = cropped_pixel_bounds.p_max.x - cropped_pixel_bounds.p_min.x;
		int offset = (p.x - cropped_pixel_bounds.p_min.x) + (p.y - cropped_pixel_bounds.p_min.y) * width;
		return pixels[offset];
	}
};

class FilmTile {
public:
	FilmTile(const Bounds2i& pb, const Vector2f& fr, const float* ft, int table_size);

	void addSample(const Point2f& p_film, const Spectrum& L, float sample_wt = 1.0f);

	FilmTilePixel& getPixel(const Point2i& p);

	Bounds2i getPixelBounds() const { return pixel_bounds; }

private:
	const Bounds2i pixel_bounds;
	const Vector2f filter_radius, inv_filter_radius;
	const float* filter_table;
	const int filter_table_size;
	std::vector<FilmTilePixel> pixels;
};

} // namespace civet

#endif // CIVET_FILM_H
