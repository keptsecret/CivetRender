#include <core/filter.h>

namespace civet {

float BoxFilter::evaluate(const Point2f& p) const {
	return 1;
}

float TriangleFilter::evaluate(const Point2f& p) const {
	return std::max(0.0f, radius.x - std::abs(p.x)) * std::max(0.0f, radius.y - std::abs(p.y));
}

float GaussianFilter::evaluate(const Point2f& p) const {
	return gaussian(p.x, exp_x) * gaussian(p.y, exp_y);
}

float GaussianFilter::gaussian(float d, float exp_v) const {
	return std::max(0.0f, std::exp(-alpha * d * d) - exp_v);
}

float MitchellFilter::evaluate(const Point2f& p) const {
	return mitchell1D(p.x * inv_radius.x) * mitchell1D(p.y * inv_radius.y);
}

float MitchellFilter::mitchell1D(float x) const {
	x = std::abs(2 * x);
	if (x > 1) {
		return ((-B - 6 * C) * x * x * x + (6 * B + 30 * C) * x * x +
					   (-12 * B - 48 * C) * x + (8 * B + 24 * C)) *
				(1.f / 6.f);
	} else {
		return ((12 - 9 * B - 6 * C) * x * x * x +
					   (-18 + 12 * B + 6 * C) * x * x +
					   (6 - 2 * B)) *
				(1.f / 6.f);
	}
}

float SincFilter::evaluate(const Point2f& p) const {
	return windowedSinc(p.x, radius.x) * windowedSinc(p.y, radius.y);
}

} // namespace civet