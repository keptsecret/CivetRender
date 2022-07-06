#ifndef CIVET_H
#define CIVET_H

#include <stdint.h>

#ifndef CIVET_L1_CACHE_LINE_SIZE
#define CIVET_L1_CACHE_LINE_SIZE 64
#endif

#if defined(__CUDA_ARCH__)
#define IS_GPU_CODE
#endif

#if defined(CIVET_BUILD_GPU) && defined(__CUDACC__)
#include <cuda.h>
#include <cuda_runtime_api.h>
#define CIVET_CPU_GPU __host__ __device__
#define CIVET_GPU __device__
#if defined(IS_GPU_CODE)
#define CIVET_CONST __device__ const
#else
#define CIVET_CONST const
#endif
#else
#define CIVET_CONST const
#define CIVET_CPU_GPU
#define CIVET_GPU
#endif

#include <float.h>
#include <limits.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <new>
#include <ostream>
#include <string>
#include <vector>
#ifdef CIVET_HAVE_MALLOC_H
#include <malloc.h>  // for _alloca, memalign
#endif
#ifdef CIVET_HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <thread>
#include <type_traits>
#include <typeinfo>
#include <utility>

#define ALLOCA(TYPE, COUNT) (TYPE*)alloca((COUNT) * sizeof(TYPE))

namespace civet {

// Global forward declarations
template <typename T>
class Vector2;
template <typename T>
class Vector3;
template <typename T>
class Point3;
template <typename T>
class Point2;
template <typename T>
class Normal3;
using Point2f = Point2<float>;
using Point2i = Point2<int>;
using Point3f = Point3<float>;
using Vector2f = Vector2<float>;
using Vector2i = Vector2<int>;
using Vector3f = Vector3<float>;
using Normal3f = Normal3<float>;

template <typename T>
class Bounds2;
using Bounds2f = Bounds2<float>;
using Bounds2i = Bounds2<int>;
template <typename T>
class Bounds3;
using Bounds3f = Bounds3<float>;
using Bounds3i = Bounds3<int>;
class Transform;

class Ray;
class RayDifferential;
struct Interaction;
class SurfaceInteraction;

class Shape;
class Primitive;

template <typename T>
class Texture;
class BSDF;
class BSSRDF;

class Medium;
class MediumInteraction;
class MediumInterface;

// Global constants
#ifdef _MSC_VER
#define MaxFloat std::numeric_limits<float>::max()
#define Infinity std::numeric_limits<float>::infinity()
#else
static CIVET_CONSTEXPR float MaxFloat = std::numeric_limits<float>::max();
static CIVET_CONSTEXPR float Infinity = std::numeric_limits<float>::infinity();
#endif
#ifdef _MSC_VER
#define MachineEpsilon (std::numeric_limits<float>::epsilon() * 0.5)
#else
static constexpr float MachineEpsilon =
		std::numeric_limits<float>::epsilon() * 0.5;
#endif
static CIVET_CONSTEXPR float ShadowEpsilon = 0.0001f;
static CIVET_CONSTEXPR float Pi = 3.14159265358979323846;
static CIVET_CONSTEXPR float InvPi = 0.31830988618379067154;
static CIVET_CONSTEXPR float Inv2Pi = 0.15915494309189533577;
static CIVET_CONSTEXPR float Inv4Pi = 0.07957747154594766788;
static CIVET_CONSTEXPR float PiOver2 = 1.57079632679489661923;
static CIVET_CONSTEXPR float PiOver4 = 0.78539816339744830961;
static CIVET_CONSTEXPR float Sqrt2 = 1.41421356237309504880;

// Global inline functions
template <typename T>
CIVET_CPU_GPU inline void swap(T& a, T& b) {
	T tmp = std::move(a);
	a = std::move(b);
	b = std::move(tmp);
}

CIVET_CPU_GPU
inline float lerp(float t, float v1, float v2) {
	return (1 - t) * v1 + t * v2;
}

CIVET_CPU_GPU
inline bool quadratic(const float a, const float b, const float c, float& t0, float& t1) {
	// find discriminant
	double discriminant = double(b) * double(b) - 4 * double(a) * double(c);
	if (discriminant < 0) {
		return false;
	}
	double rt_discrim = std::sqrt(discriminant);

	// compute t values
	double q;
	if (b < 0) {
		q = -0.5 * (b - rt_discrim);
	} else {
		q = -0.5 * (b + rt_discrim);
	}
	t0 = q / a;
	t1 = c / q;
	if (t0 > t1) {
		swap(t0, t1);
	}
	return true;
}

CIVET_CPU_GPU
inline float gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

template <typename T, typename U, typename V>
CIVET_CPU_GPU inline T clamp(T val, U low, V high) {
	if (val < low) {
		return low;
	} else if (val > high) {
		return high;
	} else {
		return val;
	}
}

template <typename T>
CIVET_CPU_GPU inline T mod(T a, T b) {
	T result = a - (a / b) * b;
	return (T)((result < 0) ? result + b : result);
}

template <>
CIVET_CPU_GPU inline float mod(float a, float b) {
	return fmod(a, b);
}

CIVET_CPU_GPU
inline float radians(float deg) {
	return (Pi / 180) * deg;
}

CIVET_CPU_GPU
inline float degrees(float rad) {
	return (180 / Pi) * rad;
}

// Math functions to support both CPU and GPU

template <typename T>
inline CIVET_CPU_GPU typename std::enable_if_t<std::is_floating_point_v<T>, bool> isNaN(T v) {
#ifdef IS_GPU_CODE
	return isnan(v);
#else
	return std::isnan(v);
#endif
}

template <typename T>
inline CIVET_CPU_GPU typename std::enable_if_t<std::is_integral_v<T>, bool> isNaN(T v) {
	return false;
}

template <typename T>
inline CIVET_CPU_GPU typename std::enable_if_t<std::is_floating_point_v<T>, bool> isInf(T v) {
#ifdef IS_GPU_CODE
	return isinf(v);
#else
	return std::isinf(v);
#endif
}

template <typename T>
inline CIVET_CPU_GPU typename std::enable_if_t<std::is_integral_v<T>, bool> isInf(T v) {
	return false;
}

template <typename T>
inline CIVET_CPU_GPU CIVET_CONSTEXPR bool isPowerOf2(T v) {
	return v && !(v & (v - 1));
}

} // namespace civet

#endif // CIVET_H