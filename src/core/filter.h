#ifndef CIVET_FILTER_H
#define CIVET_FILTER_H

#include <core/civet.h>
#include <core/geometry/vecmath.h>

namespace civet {

class Filter {
public:
	Filter(const Vector2f& r) :
			radius(r), inv_radius(Vector2f(1.0f / r.x, 1.0f / r.y)) {}

	virtual float evaluate(const Point2f& p) const = 0;

	const Vector2f radius, inv_radius;
};

class BoxFilter : public Filter {
public:
	BoxFilter(const Vector2f& r) :
			Filter(r) {}

	float evaluate(const Point2f& p) const override;
};

class TriangleFilter : public Filter {
public:
	TriangleFilter(const Vector2f& r) :
			Filter(r) {}

	float evaluate(const Point2f& p) const override;
};

class GaussianFilter : public Filter {
public:
	GaussianFilter(const Vector2f& r, float alpha) :
			Filter(r), alpha(alpha), exp_x(std::exp(-alpha * r.x * r.x)), exp_y(std::exp(-alpha * r.y * r.y)) {}

	float evaluate(const Point2f& p) const override;

private:
	float gaussian(float d, float exp_v) const;

	const float alpha;
	const float exp_x, exp_y;
};

class MitchellFilter : public Filter {
public:
	MitchellFilter(const Vector2f& r, float b, float c) :
			Filter(r), B(b), C(c) {}

	float evaluate(const Point2f& p) const override;

	float mitchell1D(float x) const;

private:
	const float B, C;
};

class SincFilter : public Filter {
public:
	SincFilter(const Vector2f& r, float t) :
			Filter(r), tau(t) {}

	float evaluate(const Point2f& p) const override;

	float sinc(float x) const {
		x = std::abs(x);
		if (x < 1e-5f) {
			return 1;
		}
		return std::sin(Pi * x) / (Pi * x);
	}

	float windowedSinc(float x, float radius) const {
		x = std::abs(x);
		if (x > radius) {
			return 0;
		}
		float lanczos = sinc(x / tau);
		return sinc(x) * lanczos;
	}

private:
	const float tau;
};

} // namespace civet

#endif // CIVET_FILTER_H
