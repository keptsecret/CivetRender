#ifndef CIVET_VECMATH_H
#define CIVET_VECMATH_H

#include <core/civet.h>

namespace civet {

/* -------------------------------------------------------------------------------
 * Vectors
 */

template <typename T>
class Vector2 {
public:
	Vector2() :
			x(0), y(0) {}

	CIVET_CPU_GPU
	Vector2(T _x, T _y) :
			x(_x), y(_y) {
		if (hasNaNs()) {
			printf("ERROR::Vector2: NaN value assigned on initialization.\n");
		}
	}

	CIVET_CPU_GPU explicit Vector2(const Point2<T>& p);
	CIVET_CPU_GPU explicit Vector2(const Point3<T>& p);

	CIVET_CPU_GPU
	Vector2<T> operator+(const Vector2<T>& v) const {
		return Vector2(x + v.x, y + v.y);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator+=(const Vector2<T>& v) const {
		x += v.x;
		y += v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator-(const Vector2<T>& v) const {
		return Vector2(x - v.x, y - v.y);
	}

	CIVET_CPU_GPU
	Vector2<T> operator-() const {
		return Vector2(-x, -y);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator-=(const Vector2<T>& v) const {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator*(T s) const {
		return Vector2(s * x, s * y);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator*=(T s) {
		x *= s;
		y *= s;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator/(T s) const {
		if (s == 0.0f) {
			printf("ERROR::Vector2: Divide by zero.\n");
		}
		float inv = 1.0f / s;
		return Vector2(x * inv, y * inv);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator/=(T s) const {
		if (s == 0.0f) {
			printf("ERROR::Vector2: Divide by zero.\n");
		}
		float inv = 1.0f / s;
		x *= inv;
		y *= inv;
		return *this;
	}

	CIVET_CPU_GPU
	bool operator==(const Vector2<T>& v) const {
		return x == v.x && y == v.y;
	}

	CIVET_CPU_GPU
	bool operator!=(const Vector2<T>& v) const {
		return x != v.x || y != v.y;
	}

	CIVET_CPU_GPU
	T operator[](int i) const {
		if (i < 0 || i > 1) {
			printf("ERROR::Vector2: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 1) {
			printf("ERROR::Vector2: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	float lengthSquared() const { return x * x + y * y; }
	CIVET_CPU_GPU
	float length() const { return std::sqrt(lengthSquared()); }

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return civet::isNaN(x) || civet::isNaN(y);
	}

	CIVET_CPU_GPU
	friend std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
		os << "[ " << v.x << ", " << v.y << " ]";
		return os;
	}

	T x, y;
};

template <typename T>
class Vector3 {
public:
	Vector3() :
			x(0), y(0), z(0) {}

	CIVET_CPU_GPU
	Vector3(T _x, T _y, T _z) :
			x(_x), y(_y), z(_z) {
		if (hasNaNs()) {
			printf("ERROR::Vector3: NaN value assigned on initialization.\n");
		}
	}

	CIVET_CPU_GPU explicit Vector3(const Point3<T>& p);
	CIVET_CPU_GPU explicit Vector3(const Normal3<T>& n);

	CIVET_CPU_GPU
	Vector3<T> operator+(const Vector3<T>& v) const {
		return Vector3(x + v.x, y + v.y, z + v.z);
	}

	CIVET_CPU_GPU
	Vector3<T>& operator+=(const Vector3<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	CIVET_CPU_GPU
	Vector3<T> operator-(const Vector3<T>& v) const {
		return Vector3(x - v.x, y - v.y, z - v.z);
	}

	CIVET_CPU_GPU
	Vector3<T> operator-() const {
		return Vector3(-x, -y, -z);
	}

	CIVET_CPU_GPU
	Vector3<T>& operator-=(const Vector3<T>& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	CIVET_CPU_GPU
	Vector3<T> operator*(T s) const {
		return Vector3(s * x, s * y, s * z);
	}

	CIVET_CPU_GPU
	Vector3<T>& operator*=(T s) {
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}

	CIVET_CPU_GPU
	Vector3<T> operator/(T s) const {
		if (s == 0.0f) {
			printf("ERROR::Vector3: Divide by zero.\n");
		}
		float inv = 1.0f / s;
		return Vector3(x * inv, y * inv, z * inv);
	}

	CIVET_CPU_GPU
	Vector3<T>& operator/=(T s) {
		if (s == 0.0f) {
			printf("ERROR::Vector3: Divide by zero.\n");
		}
		float inv = 1.0f / s;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}

	CIVET_CPU_GPU
	bool operator==(const Vector3<T>& v) const {
		return x == v.x && y == v.y && z == v.z;
	}

	CIVET_CPU_GPU
	bool operator!=(const Vector3<T>& v) const {
		return x != v.x || y != v.y || z != v.z;
	}

	CIVET_CPU_GPU
	T operator[](int i) const {
		if (i < 0 || i > 2) {
			printf("ERROR::Vector3: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 2) {
			printf("ERROR::Vector3: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	CIVET_CPU_GPU
	float lengthSquared() const { return x * x + y * y + z * z; }
	CIVET_CPU_GPU
	float length() const { return std::sqrt(lengthSquared()); }

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return civet::isNaN(x) || civet::isNaN(y) || civet::isNaN(z);
	}

	T x, y, z;
};

/* -------------------------------------------------------------------------------
 * Point
 */

template <typename T>
class Point3 {
public:
	Point3<T>() :
			x(0), y(0), z(0) {}

	CIVET_CPU_GPU
	Point3<T>(T _x, T _y, T _z) :
			x(_x), y(_y), z(_z) {
		if (hasNaNs()) {
			printf("ERROR::Point3: NaN value assigned on initialization.\n");
		}
	}

	template <typename U>
	explicit Point3(const Point3<U>& p) :
			x((T)p.x), y((T)p.y), z((T)p.z) {
		printf("ERROR::Point3: NaN value assigned on initialization.\n");
	}

	template <typename U>
	CIVET_CPU_GPU explicit operator Vector3<U>() const {
		return Vector3<U>(x, y, z);
	}

	CIVET_CPU_GPU
	Point3<T> operator+(const Vector3<T>& v) const {
		return Point3<T>(x + v.x, y + v.y, z + v.z);
	}

	CIVET_CPU_GPU
	Point3<T>& operator+=(const Vector3<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	CIVET_CPU_GPU
	Point3<T> operator+(const Point3<T>& v) const {
		return Point3<T>(x + v.x, y + v.y, z + v.z);
	}

	CIVET_CPU_GPU
	Point3<T>& operator+=(const Point3<T>& v) {
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}

	CIVET_CPU_GPU
	Vector3<T> operator-(const Point3<T>& p) const {
		return Vector3<T>(x - p.x, y - p.y, z - p.z);
	}

	CIVET_CPU_GPU
	Point3<T> operator-(const Vector3<T>& v) const {
		return Point3<T>(x - v.x, y - v.y, z - v.z);
	}

	CIVET_CPU_GPU
	Point3<T>& operator-=(const Vector3<T>& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	template <typename U>
	CIVET_CPU_GPU Point3<T> operator*(U f) const {
		return Point3<T>(f * x, f * y, f * z);
	}

	template <typename U>
	CIVET_CPU_GPU Point3<T>& operator*=(U f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	template <typename U>
	CIVET_CPU_GPU Point3<T> operator/(U f) const {
		if (f == 0.0f) {
			printf("ERROR::Point3: Divide by zero.\n");
		}
		float inv = 1.0f / f;
		return Point3<T>(inv * x, inv * y, inv * z);
	}

	template <typename U>
	CIVET_CPU_GPU Point3<T>& operator/=(U f) {
		if (f == 0.0f) {
			printf("ERROR::Point3: Divide by zero.\n");
		}
		float inv = 1.0f / f;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}

	CIVET_CPU_GPU
	T operator[](int i) const {
		if (i < 0 || i > 2) {
			printf("ERROR::Point3: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 2) {
			printf("ERROR::Point3: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	CIVET_CPU_GPU
	bool operator==(const Point3<T>& p) const {
		return x == p.x && y == p.y && z == p.z;
	}

	CIVET_CPU_GPU
	bool operator!=(const Point3<T>& p) const {
		return x != p.x || y != p.y || z != p.z;
	}

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return civet::isNaN(x) || civet::isNaN(y) || civet::isNaN(z);
	}

	T x, y, z;
};

template <typename T>
class Point2 {
public:
	Point2<T>() :
			x(0), y(0) {}

	CIVET_CPU_GPU
	Point2<T>(T _x, T _y) :
			x(_x), y(_y) {
		if (hasNaNs()) {
			printf("ERROR::Point2: NaN value assigned on initialization.\n");
		}
	}

	CIVET_CPU_GPU
	explicit Point2(const Point3<T>& p) :
			x(p.x), y(p.y) {
		if (hasNaNs()) {
			printf("ERROR::Point2: NaN value assigned on initialization.\n");
		}
	}

	template <typename U>
	explicit Point2(const Point2<U>& p) {
		x = (T)p.x;
		y = (T)p.y;
	}

	template <typename U>
	explicit Point2(const Vector2<U>& p) {
		x = (T)p.x;
		y = (T)p.y;
	}

	template <typename U>
	explicit operator Vector2<U>() const {
		return Vector2<U>(x, y);
	}

	CIVET_CPU_GPU
	Point2<T> operator+(const Vector2<T>& v) const {
		return Point2<T>(x + v.x, y + v.y);
	}

	CIVET_CPU_GPU
	Point2<T>& operator+=(const Vector2<T>& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Point2<T> operator+(const Point2<T>& v) const {
		return Point2<T>(x + v.x, y + v.y);
	}

	CIVET_CPU_GPU
	Point2<T>& operator+=(const Point2<T>& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator-(const Point2<T>& p) const {
		return Vector2<T>(x - p.x, y - p.y);
	}

	CIVET_CPU_GPU
	Point2<T> operator-(const Vector2<T>& v) const {
		return Point2<T>(x - v.x, y - v.y);
	}

	CIVET_CPU_GPU
	Point2<T>& operator-=(const Vector2<T>& v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	template <typename U>
	CIVET_CPU_GPU Point2<T> operator*(U f) const {
		return Point2<T>(f * x, f * y);
	}

	template <typename U>
	CIVET_CPU_GPU Point2<T>& operator*=(U f) {
		x *= f;
		y *= f;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator/(T s) const {
		if (s == 0.0f) {
			printf("ERROR::Point2: Divide by zero.\n");
		}
		float inv = 1 / s;
		return Vector2(x * inv, y * inv);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator/=(T s) const {
		if (s == 0.0f) {
			printf("ERROR::Point2: Divide by zero.\n");
		}
		float inv = 1 / s;
		x *= inv;
		y *= inv;
		return *this;
	}

	CIVET_CPU_GPU
	bool operator==(const Point2<T>& v) const {
		return x == v.x && y == v.y;
	}

	CIVET_CPU_GPU
	bool operator!=(const Point2<T>& v) const {
		return x != v.x || y != v.y;
	}

	CIVET_CPU_GPU
	T operator[](int i) const {
		if (i < 0 || i > 1) {
			printf("ERROR::Point2: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 1) {
			printf("ERROR::Point2: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return civet::isNaN(x) || civet::isNaN(y);
	}

	T x, y;
};

/* -------------------------------------------------------------------------------
 * Normal
 */

template <typename T>
class Normal3 {
public:
	Normal3<T>() :
			x(0), y(0), z(0) {}

	CIVET_CPU_GPU
	Normal3<T>(T _x, T _y, T _z) :
			x(_x), y(_y), z(_z) {
		if (hasNaNs()) {
			printf("ERROR::Normal3: NaN value assigned on initialization.\n");
		}
	}

	CIVET_CPU_GPU
	explicit Normal3<T>(const Vector3<T>& v) :
			x(v.x), y(v.y), z(v.z) {
		if (hasNaNs()) {
			printf("ERROR::Normal3: NaN value assigned on initialization.\n");
		}
	}

	CIVET_CPU_GPU
	Normal3<T> operator+(const Normal3<T>& n) const {
		return Normal3<T>(x + n.x, y + n.y, z + n.z);
	}

	CIVET_CPU_GPU
	Normal3<T>& operator+=(const Normal3<T>& n) {
		x += n.x;
		y += n.y;
		z += n.z;
		return *this;
	}

	CIVET_CPU_GPU
	Normal3<T> operator-(const Normal3<T>& n) const {
		return Normal3<T>(x - n.x, y - n.y, z - n.z);
	}

	CIVET_CPU_GPU
	Normal3<T> operator-() const { return Normal3(-x, -y, -z); }

	CIVET_CPU_GPU
	Normal3<T>& operator-=(const Normal3<T>& n) {
		x -= n.x;
		y -= n.y;
		z -= n.z;
		return *this;
	}

	template <typename U>
	CIVET_CPU_GPU Normal3<T> operator*(U f) const {
		return Normal3<T>(f * x, f * y, f * z);
	}

	template <typename U>
	CIVET_CPU_GPU Normal3<T>& operator*=(U f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	template <typename U>
	CIVET_CPU_GPU Normal3<T> operator/(U f) const {
		float inv = 1.0f / f;
		return Normal3<T>(x * inv, y * inv, z * inv);
	}

	template <typename U>
	CIVET_CPU_GPU Normal3<T>& operator/=(U f) {
		CHECK_NE(f, 0);
		float inv = 1.0f / f;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}

	CIVET_CPU_GPU
	bool operator==(const Normal3<T>& n) const {
		return x == n.x && y == n.y && z == n.z;
	}

	CIVET_CPU_GPU
	bool operator!=(const Normal3<T>& n) const {
		return x != n.x || y != n.y || z != n.z;
	}

	CIVET_CPU_GPU
	T operator[](int i) const {
		if (i < 0 || i > 2) {
			printf("ERROR::Normal3: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 2) {
			printf("ERROR::Normal3: Tried to access element at index %d.\n", i);
		}
		if (i == 0) {
			return x;
		}
		if (i == 1) {
			return y;
		}
		return z;
	}

	CIVET_CPU_GPU
	float lengthSquared() const { return x * x + y * y + z * z; }
	CIVET_CPU_GPU
	float length() const { return std::sqrt(lengthSquared()); }

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return civet::isNaN(x) || civet::isNaN(y) || civet::isNaN(z);
	}

	T x, y, z;
};

/* -------------------------------------------------------------------------------
 * Bounds (AABBs)
 */

template <typename T>
class Bounds2 {
public:
	Bounds2() {
		T min_num = std::numeric_limits<T>::lowest();
		T max_num = std::numeric_limits<T>::max();
		p_min = Point2<T>(max_num, max_num);
		p_max = Point2<T>(min_num, min_num);
	}

	CIVET_CPU_GPU
	explicit Bounds2(const Point2<T>& p) :
			p_min(p), p_max(p) {}

	template <typename U>
	explicit operator Bounds2<U>() const {
		return Bounds2<U>((Point2<U>)p_min, (Point2<U>)p_max);
	}

	CIVET_CPU_GPU
	Bounds2(const Point2<T>& p1, const Point2<T>& p2) :
			p_min(std::min(p1.x, p2.x), std::min(p1.y, p2.y)), p_max(std::max(p1.x, p2.x), std::max(p1.y, p2.y)) {}

	CIVET_CPU_GPU
	inline const Point2<T>& operator[](int i) const {
		return (i == 0) ? p_min : p_max;
	}

	CIVET_CPU_GPU
	inline Point2<T>& operator[](int i) {
		return (i == 0) ? p_min : p_max;
	}

	bool operator==(const Bounds2<T>& b) const {
		return b.p_min == p_min && b.p_max == p_max;
	}

	bool operator!=(const Bounds2<T>& b) const {
		return b.p_min != p_min || b.p_max != p_max;
	}

	CIVET_CPU_GPU
	Vector2<T> diagonal() const { return p_max - p_min; }

	CIVET_CPU_GPU
	T area() const {
		Vector2<T> d = diagonal();
		return d.x * d.y;
	}

	CIVET_CPU_GPU
	int maximumAxis() const {
		Vector2<T> d = diagonal();
		if (d.x > d.y) {
			return 0;
		} else {
			return 1;
		}
	}

	CIVET_CPU_GPU
	Point2<T> lerp(const Point2f& t) const {
		return Point2<T>(civet::lerp(t.x, p_min.x, p_max.x), civet::lerp(t.y, p_min.y, p_max.y));
	}

	CIVET_CPU_GPU
	Vector2<T> offset(const Point2<T>& p) const {
		Vector2<T> o = p - p_min;
		if (p_max.x > p_min.x) {
			o.x /= p_max.x - p_min.x;
		}
		if (p_max.y > p_min.y) {
			o.y /= p_max.y - p_min.y;
		}
		return o;
	}

	CIVET_CPU_GPU
	void boundingSphere(Point2<T>* center, float* radius) const {
		*center = (p_min + p_max) / 2;
		*radius = bInside(*center, *this) ? distance(*center, p_max) : 0;
	}

	CIVET_CPU_GPU
	friend std::ostream& operator<<(std::ostream& os, const Bounds2<T>& b) {
		os << "[ " << b.pMin << " - " << b.pMax << " ]";
		return os;
	}

	Point2<T> p_min, p_max;
};

template <typename T>
class Bounds3 {
public:
	Bounds3() {
		T min_num = std::numeric_limits<T>::lowest();
		T max_num = std::numeric_limits<T>::max();
		p_min = Point3<T>(max_num, max_num, max_num);
		p_max = Point3<T>(min_num, min_num, min_num);
	}

	CIVET_CPU_GPU
	explicit Bounds3(const Point3<T>& p) :
			p_min(p), p_max(p) {}

	CIVET_CPU_GPU
	Bounds3(const Point3<T>& p1, const Point3<T>& p2) :
			p_min(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z)), p_max(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z)) {}

	CIVET_CPU_GPU
	inline const Point3<T>& operator[](int i) const {
		return (i == 0) ? p_min : p_max;
	}

	CIVET_CPU_GPU
	inline Point3<T>& operator[](int i) {
		return (i == 0) ? p_min : p_max;
	}

	bool operator==(const Bounds3<T>& b) const {
		return b.p_min == p_min && b.p_max == p_max;
	}

	bool operator!=(const Bounds3<T>& b) const {
		return b.p_min != p_min || b.p_max != p_max;
	}

	CIVET_CPU_GPU
	Point3<T> corner(int n_corner) const {
		return Point3<T>((*this)[(n_corner & 1)].x,
				(*this)[(n_corner & 2) ? 1 : 0].y,
				(*this)[(n_corner & 4) ? 1 : 0].z);
	}

	CIVET_CPU_GPU
	Vector3<T> diagonal() const { return p_max - p_min; }

	CIVET_CPU_GPU
	T surfaceArea() const {
		Vector3<T> d = diagonal();
		return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
	}

	CIVET_CPU_GPU
	T volume() const {
		Vector3<T> d = diagonal();
		return d.x * d.y * d.z;
	}

	CIVET_CPU_GPU
	int maximumAxis() const {
		Vector3<T> d = diagonal();
		if (d.x > d.y && d.x > d.z) {
			return 0;
		} else if (d.y > d.z) {
			return 1;
		} else {
			return 2;
		}
	}

	CIVET_CPU_GPU
	Point3<T> lerp(const Point3f& t) const {
		return Point3<T>(civet::lerp(t.x, p_min.x, p_max.x), civet::lerp(t.y, p_min.y, p_max.y), civet::lerp(t.z, p_min.z, p_max.z));
	}

	CIVET_CPU_GPU
	Vector3<T> offset(const Point3<T>& p) const {
		Vector3<T> o = p - p_min;
		if (p_max.x > p_min.x) {
			o.x /= p_max.x - p_min.x;
		}
		if (p_max.y > p_min.y) {
			o.y /= p_max.y - p_min.y;
		}
		if (p_max.z > p_min.z) {
			o.z /= p_max.z - p_min.z;
		}
		return o;
	}

	CIVET_CPU_GPU
	void boundingSphere(Point3<T>* center, float* radius) const {
		*center = (p_min + p_max) / 2;
		*radius = bInside(*center, *this) ? distance(*center, p_max) : 0;
	}

	CIVET_CPU_GPU
	bool intersectP(Point3f o, Vector3f d, float t_max = Infinity, float* hitt0 = nullptr, float* hitt1 = nullptr) const;

	CIVET_CPU_GPU
	bool intersectP(Point3f o, Vector3f d, float t_max, Vector3f& inv_dir, const int dir_is_neg[3]) const;

	CIVET_CPU_GPU
	friend std::ostream& operator<<(std::ostream& os, const Bounds3<T>& b) {
		os << "[ " << b.pMin << " - " << b.pMax << " ]";
		return os;
	}

	Point3<T> p_min, p_max;
};

class Bounds2iIterator : public std::forward_iterator_tag {
public:
	CIVET_CPU_GPU
	Bounds2iIterator(const Bounds2i& b, const Point2i& pt) :
			p(pt), bounds(&b) {}

	CIVET_CPU_GPU
	Bounds2iIterator operator++() {
		advance();
		return *this;
	}

	CIVET_CPU_GPU
	Bounds2iIterator operator++(int) {
		Bounds2iIterator old = *this;
		advance();
		return old;
	}

	CIVET_CPU_GPU
	bool operator==(const Bounds2iIterator& bi) const {
		return p == bi.p && bounds == bi.bounds;
	}

	CIVET_CPU_GPU
	bool operator!=(const Bounds2iIterator& bi) const {
		return p != bi.p || bounds != bi.bounds;
	}

	CIVET_CPU_GPU
	Point2i operator*() const { return p; }

private:
	CIVET_CPU_GPU
	void advance() {
		++p.x;
		if (p.x == bounds->p_max.x) {
			p.x = bounds->p_min.x;
			++p.y;
		}
	}

	Point2i p;
	const Bounds2i* bounds;
};

inline Bounds2iIterator begin(const Bounds2i &b) {
	return Bounds2iIterator(b, b.p_min);
}

inline Bounds2iIterator end(const Bounds2i &b) {
	// Normally, the ending point is at the minimum x value and one past
	// the last valid y value.
	Point2i pEnd(b.p_min.x, b.p_max.y);
	// However, if the bounds are degenerate, override the end point to
	// equal the start point so that any attempt to iterate over the bounds
	// exits out immediately.
	if (b.p_min.x >= b.p_max.x || b.p_min.y >= b.p_max.y)	{
		pEnd = b.p_min;
	}
	return Bounds2iIterator(b, pEnd);
}

/*---------------------------------------------------------------------------------*/
/*
 * Inline geometry functions
 */

// Vector2
template <typename T>
CIVET_CPU_GPU Vector2<T>::Vector2(const Point2<T>& p) :
		x(p.x), y(p.y) {
	if (hasNaNs()) {
		printf("ERROR::Vector2: NaN value assigned on initialization.\n");
	}
}

template <typename T>
CIVET_CPU_GPU Vector2<T>::Vector2(const Point3<T>& p) :
		x(p.x), y(p.y) {
	if (hasNaNs()) {
		printf("ERROR::Vector2: NaN value assigned on initialization.\n");
	}
}

template <typename T>
CIVET_CPU_GPU inline Vector2<T> operator*(T s, const Vector2<T>& v) {
	return v * s;
}

template <typename T>
CIVET_CPU_GPU inline Vector2<T> abs(const Vector2<T>& v) {
	return Vector2<T>(std::abs(v.x), std::abs(v.y));
}

template <typename T>
CIVET_CPU_GPU inline T dot(const Vector2<T>& v1, const Vector2<T>& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

template <typename T>
CIVET_CPU_GPU inline T absDot(const Vector2<T>& v1, const Vector2<T>& v2) {
	return std::abs(dot(v1, v2));
}

template <typename T>
CIVET_CPU_GPU inline Vector2<T> normalize(const Vector2<T>& v) {
	return v / v.length();
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Vector2<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

// Vector3
template <typename T>
CIVET_CPU_GPU Vector3<T>::Vector3(const Point3<T>& p) :
		x(p.x), y(p.y), z(p.z) {
	if (hasNaNs()) {
		printf("ERROR::Vector3: NaN value assigned on initialization.\n");
	}
}

template <typename T>
CIVET_CPU_GPU Vector3<T>::Vector3(const Normal3<T>& n) :
		x(n.x), y(n.y), z(n.z) {
	if (hasNaNs()) {
		printf("ERROR::Vector3: NaN value assigned on initialization.\n");
	}
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> operator*(T s, const Vector3<T>& v) {
	return v * s;
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> abs(const Vector3<T>& v) {
	return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template <typename T>
CIVET_CPU_GPU inline T dot(const Vector3<T>& v1, const Vector3<T>& v2) {
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
CIVET_CPU_GPU inline T absDot(const Vector3<T>& v1, const Vector3<T>& v2) {
	return std::abs(dot(v1, v2));
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> cross(const Vector3<T>& v1, const Vector3<T>& v2) {
	double v1x = v1.x, v1y = v1.y, v1z = v1.z;
	double v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector3<T>(v1y * v2z - v1z * v2y,
			v1z * v2x - v1x * v2z,
			v1x * v2y - v1y * v2z);
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> cross(const Vector3<T>& v, const Normal3<T>& n) {
	double v1x = v.x, v1y = v.y, v1z = v.z;
	double v2x = n.x, v2y = n.y, v2z = n.z;
	return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
			(v1x * v2y) - (v1y * v2x));
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> cross(const Normal3<T>& n, const Vector3<T>& v) {
	double v1x = n.x, v1y = n.y, v1z = n.z;
	double v2x = v.x, v2y = v.y, v2z = v.z;
	return Vector3<T>((v1y * v2z) - (v1z * v2y), (v1z * v2x) - (v1x * v2z),
			(v1x * v2y) - (v1y * v2x));
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> normalize(const Vector3<T>& v) {
	return v / v.length();
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

template <typename T>
CIVET_CPU_GPU inline T minComponent(const Vector3<T>& v) {
	return std::min(v.x, std::min(v.y, v.z));
}

template <typename T>
CIVET_CPU_GPU inline T maxComponent(const Vector3<T>& v) {
	return std::max(v.x, std::max(v.y, v.z));
}

template <typename T>
CIVET_CPU_GPU inline int maxDimension(const Vector3<T>& v) {
	return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> min(const Vector3<T>& v1, const Vector3<T>& v2) {
	return Vector3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> max(const Vector3<T>& v1, const Vector3<T>& v2) {
	return Vector3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
}

template <typename T>
CIVET_CPU_GPU inline Vector3<T> permute(const Vector3<T>& v, int x, int y, int z) {
	return Vector3<T>(v[x], v[y], v[z]);
}

template <typename T>
CIVET_CPU_GPU inline void coordinateSystem(const Vector3<T>& v1, Vector3<T>* v2, Vector3<T>* v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		*v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	} else {
		*v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	}
	*v3 = cross(v1, *v2);
}

// Point2
template <typename T>
CIVET_CPU_GPU inline float distance(const Point2<T>& p1, const Point2<T>& p2) {
	return (p1 - p2).length();
}

template <typename T>
CIVET_CPU_GPU inline float distanceSquared(const Point2<T>& p1, const Point2<T>& p2) {
	return (p1 - p2).lengthSquared();
}

template <typename T, typename U>
CIVET_CPU_GPU inline Point2<T> operator*(U f, const Point2<T>& p) {
	return p * f;
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> lerp(float t, const Point2<T>& p1, const Point2<T>& p2) {
	return (1 - t) * p1 + t * p2;
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> min(const Point2<T>& v1, const Point2<T>& v2) {
	return Point2<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y));
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> max(const Point2<T>& v1, const Point2<T>& v2) {
	return Point2<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y));
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> floor(const Point2<T>& p) {
	return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> ceil(const Point2<T>& p) {
	return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> abs(const Point2<T>& p) {
	return Point2<T>(std::abs(p.x), std::abs(p.y));
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Point2<T>& v) {
	os << "[ " << v.x << ", " << v.y << " ]";
	return os;
}

// Point3
template <typename T>
CIVET_CPU_GPU inline float distance(const Point3<T>& p1, const Point3<T>& p2) {
	return (p1 - p2).length();
}

template <typename T>
CIVET_CPU_GPU inline float distanceSquared(const Point3<T>& p1, const Point3<T>& p2) {
	return (p1 - p2).lengthSquared();
}

template <typename T, typename U>
CIVET_CPU_GPU inline Point3<T> operator*(U f, const Point3<T>& p) {
	return p * f;
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> lerp(float t, const Point3<T>& p1, const Point3<T>& p2) {
	return (1 - t) * p1 + t * p2;
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> min(const Point3<T>& v1, const Point3<T>& v2) {
	return Point3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> max(const Point3<T>& v1, const Point3<T>& v2) {
	return Point3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> floor(const Point3<T>& p) {
	return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> ceil(const Point3<T>& p) {
	return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> abs(const Point3<T>& p) {
	return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> permute(const Point3<T>& p, int x, int y, int z) {
	return Point3<T>(p[x], p[y], p[z]);
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Point3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

// Normal
template <typename T, typename U>
CIVET_CPU_GPU inline Normal3<T> operator*(U f, const Normal3<T>& n) {
	return Normal3<T>(f * n.x, f * n.y, f * n.z);
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> normalize(const Normal3<T>& n) {
	return n / n.length();
}

template <typename T>
inline T dot(const Normal3<T>& n1, const Vector3<T>& v2) {
	return n1.x * v2.x + n1.y * v2.y + n1.z * v2.z;
}

template <typename T>
CIVET_CPU_GPU inline T dot(const Vector3<T>& v1, const Normal3<T>& n2) {
	return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

template <typename T>
CIVET_CPU_GPU inline T dot(const Normal3<T>& n1, const Normal3<T>& n2) {
	return n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
}

template <typename T>
CIVET_CPU_GPU inline T absDot(const Normal3<T>& n1, const Vector3<T>& v2) {
	return std::abs(n1.x * v2.x + n1.y * v2.y + n1.z * v2.z);
}

template <typename T>
CIVET_CPU_GPU inline T absDot(const Vector3<T>& v1, const Normal3<T>& n2) {
	return std::abs(v1.x * n2.x + v1.y * n2.y + v1.z * n2.z);
}

template <typename T>
CIVET_CPU_GPU inline T absDot(const Normal3<T>& n1, const Normal3<T>& n2) {
	return std::abs(n1.x * n2.x + n1.y * n2.y + n1.z * n2.z);
}

template <typename T>
CIVET_CPU_GPU inline T abs(const Normal3<T>& n) {
	return Normal3<T>(std::abs(n.x), std::abs(n.y), std::abs(n.z));
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> faceforward(const Normal3<T>& n, const Vector3<T>& v) {
	return dot(n, v) < 0.0f ? -n : n;
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> faceforward(const Normal3<T>& n, const Normal3<T>& n2) {
	return dot(n, n2) < 0.0f ? -n : n;
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> faceforward(const Vector3<T>& v, const Vector3<T>& v2) {
	return dot(v, v2) < 0.0f ? -v : v;
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> faceforward(const Vector3<T>& v, const Normal3<T>& n) {
	return dot(v, n) < 0.0f ? -v : v;
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Normal3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
}

// Bounds2
template <typename T>
CIVET_CPU_GPU inline Bounds2<T> bUnion(const Bounds2<T>& b, const Point2<T>& p) {
	return Bounds2<T>(min(b.p_min, p), max(b.p_max, p));
}

template <typename T>
CIVET_CPU_GPU inline Bounds2<T> bUnion(const Bounds2<T>& b1, const Bounds2<T>& b2) {
	return Bounds2<T>(min(b1.p_min, b2.p_min), max(b1.p_max, b2.p_max));
}

template <typename T>
CIVET_CPU_GPU inline Bounds2<T> bIntersect(const Bounds2<T>& b1, const Bounds2<T>& b2) {
	Bounds2<T> b;
	b.p_min = max(b1.p_min, b2.p_min);
	b.p_max = min(b1.p_max, b2.p_max);
	return b;
}

template <typename T>
CIVET_CPU_GPU inline bool bOverlaps(const Bounds2<T>& b1, const Bounds2<T>& b2) {
	bool x = (b1.p_max.x >= b2.p_min.x) && (b1.p_min.x <= b2.p_max.x);
	bool y = (b1.p_max.y >= b2.p_min.y) && (b1.p_min.y <= b2.p_max.y);
	return x && y;
}

template <typename T>
CIVET_CPU_GPU inline bool bInside(const Point2<T>& p, const Bounds2<T>& b) {
	return p.x >= b.p_min.x && p.x <= b.p_max.x &&
			p.y >= b.p_min.y && p.y <= b.p_max.y;
}

template <typename T>
CIVET_CPU_GPU inline bool bInsideExclusive(const Point2<T>& p, const Bounds2<T>& b) {
	return p.x >= b.p_min.x && p.x < b.p_max.x &&
			p.y >= b.p_min.y && p.y < b.p_max.y;
}

template <typename T, typename U>
CIVET_CPU_GPU inline bool bExpand(const Bounds2<T>& b, U delta) {
	return Bounds2<T>(b.p_min - Vector2<T>(delta, delta), b.p_max + Vector2<T>(delta, delta));
}

// Bounds3
template <typename T>
CIVET_CPU_GPU inline Bounds3<T> bUnion(const Bounds3<T>& b, const Point3<T>& p) {
	return Bounds3<T>(min(b.p_min, p), max(b.p_max, p));
}

template <typename T>
CIVET_CPU_GPU inline Bounds3<T> bUnion(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	return Bounds3<T>(min(b1.p_min, b2.p_min), max(b1.p_max, b2.p_max));
}

template <typename T>
CIVET_CPU_GPU inline Bounds3<T> bIntersect(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	Bounds3<T> b;
	b.p_min = max(b1.p_min, b2.p_min);
	b.p_max = min(b1.p_max, b2.p_max);
	return b;
}

template <typename T>
CIVET_CPU_GPU inline bool bOverlaps(const Bounds3<T>& b1, const Bounds3<T>& b2) {
	bool x = (b1.p_max.x >= b2.p_min.x) && (b1.p_min.x <= b2.p_max.x);
	bool y = (b1.p_max.y >= b2.p_min.y) && (b1.p_min.y <= b2.p_max.y);
	bool z = (b1.p_max.z >= b2.p_min.z) && (b1.p_min.z <= b2.p_max.z);
	return x && y && z;
}

template <typename T>
CIVET_CPU_GPU inline bool bInside(const Point3<T>& p, const Bounds3<T>& b) {
	return p.x >= b.p_min.x && p.x <= b.p_max.x &&
			p.y >= b.p_min.y && p.y <= b.p_max.y &&
			p.z >= b.p_min.z && p.z <= b.p_max.z;
}

template <typename T>
CIVET_CPU_GPU inline bool bInsideExclusive(const Point3<T>& p, const Bounds3<T>& b) {
	return p.x >= b.p_min.x && p.x < b.p_max.x &&
			p.y >= b.p_min.y && p.y < b.p_max.y &&
			p.z >= b.p_min.z && p.z < b.p_max.z;
}

template <typename T, typename U>
CIVET_CPU_GPU inline Bounds3<T> bExpand(const Bounds3<T>& b, U delta) {
	return Bounds3<T>(b.p_min - Vector3<T>(delta, delta, delta), b.p_max + Vector3<T>(delta, delta, delta));
}

template <typename T>
CIVET_CPU_GPU inline bool Bounds3<T>::intersectP(Point3f o, Vector3f d, float t_max, float* hitt0, float* hitt1) const {
	float t0 = 0, t1 = t_max;
	for (int i = 0; i < 3; i++) {
		float inv_dir = 1 / d[i];
		float t_near = (p_min[i] - o[i]) * inv_dir;
		float t_far = (p_max[i] - o[i]) * inv_dir;

		if (t_near > t_far) {
			swapElem(t_near, t_far);
		}
		t_far *= 1 + 2 * gamma(3);
		t0 = t_near > t0 ? t_near : t0;
		t1 = t_far < t1 ? t_far : t1;
		if (t0 > t1) {
			return false;
		}
	}
	if (hitt0) {
		*hitt0 = t0;
	}
	if (hitt1) {
		*hitt1 = t1;
	}
	return true;
}

template <typename T>
CIVET_CPU_GPU inline bool Bounds3<T>::intersectP(Point3f o, Vector3f d, float ray_t_max, Vector3f& inv_dir, const int dir_is_neg[3]) const {
	const Bounds3f& bounds = *this;
	float t_min = (bounds[dir_is_neg[0]].x - o.x) * inv_dir.x;
	float t_max = (bounds[1 - dir_is_neg[0]].x - o.x) * inv_dir.x;
	float ty_min = (bounds[dir_is_neg[1]].y - o.y) * inv_dir.y;
	float ty_max = (bounds[1 - dir_is_neg[1]].y - o.y) * inv_dir.y;

	t_max *= 1 + 2 * gamma(3);
	ty_max *= 1 + 2 * gamma(3);

	if (t_min > ty_max || ty_min > t_max) {
		return false;
	}
	if (ty_min > t_min) {
		t_min = ty_min;
	}
	if (ty_max < t_max) {
		t_max = ty_min;
	}

	float tz_min = (bounds[dir_is_neg[2]].z - o.z) * inv_dir.z;
	float tz_max = (bounds[1 - dir_is_neg[2]].z - o.z) * inv_dir.z;

	tz_max *= 1 + 2 * gamma(3);

	if (t_min > tz_max || tz_min > t_max) {
		return false;
	}
	if (tz_min > t_min) {
		t_min = tz_min;
	}
	if (tz_max < t_max) {
		t_max = tz_min;
	}

	return (t_min < ray_t_max) && (t_max > 0);
}

// Miscellaneous angle functions
CIVET_CPU_GPU
inline Vector3f sphericalDirection(float sin_theta, float cos_theta, float phi) {
	return Vector3f(sin_theta * cos(phi), sin_theta * sin(phi), cos_theta);
}

CIVET_CPU_GPU
inline Vector3f sphericalDirection(float sin_theta, float cos_theta, float phi, const Vector3f& x, const Vector3f& y, const Vector3f& z) {
	return sin_theta * cos(phi) * x + sin_theta * sin(phi) * y + cos_theta * z;
}

CIVET_CPU_GPU
inline float sphericalTheta(const Vector3f& v) {
	return acos(clamp(v.z, -1, 1));
}

CIVET_CPU_GPU
inline float sphericalPhi(const Vector3f& v) {
	float p = atan2(v.y, v.x);
	return (p < 0) ? (p + 2 * Pi) : p;
}

} // namespace civet

#endif // CIVET_VECMATH_H