#ifndef CIVET_GEOMETRY_H
#define CIVET_GEOMETRY_H

#include "civet.h"

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
		if (!hasNaNs()) {
			std::cerr << "ERROR::Vector2: NaN value assigned on initialization.\n";
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
	Vector2<T>& operator*=(T s) const {
		x *= s;
		y *= s;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator/(T s) const {
		if (s == 0.0f) {
			std::cerr << "ERROR::Vector2: Divide by zero.\n";
		}
		float inv = 1.0f / s;
		return Vector2(x * inv, y * inv);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator/=(T s) const {
		if (s == 0.0f) {
			std::cerr << "ERROR::Vector2: Divide by zero.\n";
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
			std::cerr << "ERROR::Vector2: Tried to access element at index " << i << ".\n";
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 1) {
			std::cerr << "ERROR::Vector2: Tried to access element at index " << i << ".\n";
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
		return std::isnan(x) || std::isnan(y);
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
		if (!hasNaNs()) {
			std::cerr << "ERROR::Vector3: NaN value assigned on initialization.\n";
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
			std::cerr << "ERROR::Vector3: Divide by zero.\n";
		}
		float inv = 1.0f / s;
		return Vector3(x * inv, y * inv, z * inv);
	}

	CIVET_CPU_GPU
	Vector3<T>& operator/=(T s) const {
		if (s == 0.0f) {
			std::cerr << "ERROR::Vector3: Divide by zero.\n";
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
			std::cerr << "ERROR::Vector3: Tried to access element at index " << i << ".\n";
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
			std::cerr << "ERROR::Vector3: Tried to access element at index " << i << ".\n";
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
		return std::isnan(x) || std::isnan(y) || std::isnan(z);
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
        if (!hasNaNs()) {
            std::cerr << "ERROR::Point3: NaN value assigned on initialization.\n";
		}
    }

    template <typename U>
    explicit Point3(const Point3<U>& p) :
            x((T)p.x), y((T)p.y), z((T)p.z) {
        std::cerr << "ERROR::Point3: NaN value assigned on initialization.\n";
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
			std::cerr << "ERROR::Point3: Divide by zero.\n";
		}
        float inv = 1.0f / f;
        return Point3<T>(inv * x, inv * y, inv * z);
    }

    template <typename U>
	CIVET_CPU_GPU Point3<T>& operator/=(U f) {
		if (f == 0.0f) {
			std::cerr << "ERROR::Point3: Divide by zero.\n";
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
			std::cerr << "ERROR::Point3: Tried to access element at index " << i << ".\n";
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
            std::cerr << "ERROR::Point3: Tried to access element at index " << i << ".\n";
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
        return std::isnan(x) || std::isnan(y) || std::isnan(z);
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
		if (!hasNaNs()) {
			std::cerr << "ERROR::Point2: NaN value assigned on initialization.\n";
		}
	}

	CIVET_CPU_GPU
	explicit Point2(const Point3<T>& p) :
			x(p.x), y(p.y) {
		if (!hasNaNs()) {
			std::cerr << "ERROR::Point2: NaN value assigned on initialization.\n";
		}
	}

	CIVET_CPU_GPU
	Point2<T> operator+(const Vector2<T>& v) const {
		return Point3<T>(x + v.x, y + v.y);
	}

	CIVET_CPU_GPU
	Point2<T>& operator+=(const Vector2<T>& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Point2<T> operator+(const Point2<T>& v) const {
		return Point3<T>(x + v.x, y + v.y);
	}

	CIVET_CPU_GPU
	Point2<T>& operator+=(const Point2<T>& v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator-(const Point2<T>& p) const {
		return Vector3<T>(x - p.x, y - p.y);
	}

	CIVET_CPU_GPU
	Point2<T> operator-(const Vector2<T>& v) const {
		return Point3<T>(x - v.x, y - v.y);
	}

	CIVET_CPU_GPU
	Point2<T>& operator-=(const Vector2<T>& v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator*(T s) const {
		return Vector2(s * x, s * y);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator*=(T s) const {
		x *= s;
		y *= s;
		return *this;
	}

	CIVET_CPU_GPU
	Vector2<T> operator/(T s) const {
		if (s == 0.0f) {
			std::cerr << "ERROR::Point2: Divide by zero.\n";
		}
		float inv = 1 / s;
		return Vector2(x * inv, y * inv);
	}

	CIVET_CPU_GPU
	Vector2<T>& operator/=(T s) const {
		if (s == 0.0f) {
			std::cerr << "ERROR::Point2: Divide by zero.\n";
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
			std::cerr << "ERROR::Point2: Tried to access element at index " << i << ".\n";
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	T& operator[](int i) {
		if (i < 0 || i > 1) {
			std::cerr << "ERROR::Point2: Tried to access element at index " << i << ".\n";
		}
		if (i == 0) {
			return x;
		}
		return y;
	}

	CIVET_CPU_GPU
	bool hasNaNs() const {
		return std::isnan(x) || std::isnan(y);
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
		if (!hasNaNs()) {
			std::cerr << "ERROR::Normal3: NaN value assigned on initialization.\n";
		}
	}

	CIVET_CPU_GPU
	explicit Normal3<T>(const Vector3<T>& v) :
	        x(v.x), y(v.y), z(v.z) {
		if (!hasNaNs()) {
			std::cerr << "ERROR::Normal3: NaN value assigned on initialization.\n";
		}
	}

	CIVET_CPU_GPU
	Normal3<T> operator+(const Normal3<T> &n) const {
		return Normal3<T>(x + n.x, y + n.y, z + n.z);
	}

	CIVET_CPU_GPU
	Normal3<T> &operator+=(const Normal3<T> &n) {
		x += n.x;
		y += n.y;
		z += n.z;
		return *this;
	}

	CIVET_CPU_GPU
	Normal3<T> operator-(const Normal3<T> &n) const {
		return Normal3<T>(x - n.x, y - n.y, z - n.z);
	}

	CIVET_CPU_GPU
	Normal3<T> operator-() const { return Normal3(-x, -y, -z); }

	CIVET_CPU_GPU
	Normal3<T> &operator-=(const Normal3<T> &n) {
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
	CIVET_CPU_GPU Normal3<T> &operator*=(U f) {
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
	CIVET_CPU_GPU Normal3<T> &operator/=(U f) {
		CHECK_NE(f, 0);
		float inv = 1.0f / f;
		x *= inv;
		y *= inv;
		z *= inv;
		return *this;
	}

	CIVET_CPU_GPU
	bool operator==(const Normal3<T> &n) const {
		return x == n.x && y == n.y && z == n.z;
	}

	CIVET_CPU_GPU
	bool operator!=(const Normal3<T> &n) const {
		return x != n.x || y != n.y || z != n.z;
	}

	CIVET_CPU_GPU
	T operator[](int i) const {
		if (i < 0 || i > 2) {
			std::cerr << "ERROR::Normal3: Tried to access element at index " << i << ".\n";
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
	T &operator[](int i) {
		if (i < 0 || i > 2) {
			std::cerr << "ERROR::Normal3: Tried to access element at index " << i << ".\n";
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
		return std::isnan(x) || std::isnan(y) || std::isnan(z);
	}

	T x, y, z;
};

// Some common typedefs
typedef Vector2<float> Vector2f;
typedef Vector2<int> Vector2i;
typedef Vector3<float> Vector3f;
typedef Vector3<int> Vector3i;

typedef Point2<float> Point2f;
typedef Point2<int> Point2i;
typedef Point3<float> Point3f;
typedef Point3<int> Point3i;

typedef Normal3<float> Normal3f;

/*---------------------------------------------------------------------------------*/
/*
 * Inline geometry functions
 */

// Vector2
template <typename T>
Vector2<T>::Vector2(const Point2<T> &p)
		: x(p.x), y(p.y) {
	if (!hasNaNs()) {
		std::cerr << "ERROR::Vector2: NaN value assigned on initialization.\n";
	}
}

template <typename T>
Vector2<T>::Vector2(const Point3<T> &p)
		: x(p.x), y(p.y) {
	if (!hasNaNs()) {
		std::cerr << "ERROR::Vector2: NaN value assigned on initialization.\n";
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
	os << "[ " << v.x << ", " << v.y << " ]";
	return os;
}

// Vector3
template <typename T>
Vector3<T>::Vector3(const Point3<T>& p) :
		x(p.x), y(p.y), z(p.z) {
	if (!hasNaNs()) {
		std::cerr << "ERROR::Vector3: NaN value assigned on initialization.\n";
	}
}

template <typename T>
Vector3<T>::Vector3(const Normal3<T>& n) :
		x(n.x), y(n.y), z(n.z) {
	if (!hasNaNs()) {
		std::cerr << "ERROR::Vector3: NaN value assigned on initialization.\n";
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
CIVET_CPU_GPU inline Vector3<T> normalize(const Vector3<T>& v) {
	return v / v.length();
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
CIVET_CPU_GPU inline Vector3<T> coordinateSystem(const Vector3<T>& v1, Vector3<T>* v2, Vector3<T>* v3) {
	if (std::abs(v1.x) > std::abs(v1.y)) {
		*v2 = Vector3<T>(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	} else {
		*v2 = Vector3<T>(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	}
	*v3 = cross(v1, *v2);
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
	os << "[ " << v.x << ", " << v.y << ", " << v.z << " ]";
	return os;
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
CIVET_CPU_GPU inline Point2<T> floor(const Point2<T> &p) {
	return Point2<T>(std::floor(p.x), std::floor(p.y));
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> ceil(const Point2<T> &p) {
	return Point2<T>(std::ceil(p.x), std::ceil(p.y));
}

template <typename T>
CIVET_CPU_GPU inline Point2<T> abs(const Point2<T> &p) {
	return Point2<T>(std::abs(p.x), std::abs(p.y));
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Point2<T>& p) {
	os << "[ " << p.x << ", " << p.y << " ]";
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
CIVET_CPU_GPU inline Point3<T> floor(const Point3<T> &p) {
    return Point3<T>(std::floor(p.x), std::floor(p.y), std::floor(p.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> ceil(const Point3<T> &p) {
    return Point3<T>(std::ceil(p.x), std::ceil(p.y), std::ceil(p.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> abs(const Point3<T> &p) {
	return Point3<T>(std::abs(p.x), std::abs(p.y), std::abs(p.z));
}

template <typename T>
CIVET_CPU_GPU inline Point3<T> permute(const Point3<T>& p, int x, int y, int z) {
	return Point3<T>(p[x], p[y], p[z]);
}

template <typename T>
CIVET_CPU_GPU inline std::ostream& operator<<(std::ostream& os, const Point3<T>& p) {
	os << "[ " << p.x << ", " << p.y << ", " << p.z << " ]";
	return os;
}

// Normal
template <typename T, typename U>
CIVET_CPU_GPU inline Normal3<T> operator*(U f, const Normal3<T>& n) {
	return Normal3<T>(f * n.x, f * n.y, f * n.z);
}

template <typename T>
CIVET_CPU_GPU inline Normal3<T> normalize(const Normal3<T>& n) {
	return n / n.Length();
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

} // namespace civet

#endif // CIVET_GEOMETRY_H