#ifndef CIVET_H
#define CIVET_H

#include <stdint.h>
#include <cstddef>

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
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <new>
#include <string>
#include <thread>
#include <type_traits>
#include <typeinfo>
#include <utility>

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
class Medium;

// Global constants
#ifdef _MSC_VER
#define MaxFloat std::numeric_limits<float>::max()
#define Infinity std::numeric_limits<float>::infinity()
#else
static constexpr float MaxFloat = std::numeric_limits<float>::max();
static constexpr float Infinity = std::numeric_limits<float>::infinity();
#endif
#ifdef _MSC_VER
#define MachineEpsilon (std::numeric_limits<float>::epsilon() * 0.5)
#else
static constexpr float MachineEpsilon =
		std::numeric_limits<float>::epsilon() * 0.5;
#endif
static constexpr float ShadowEpsilon = 0.0001f;
static constexpr float Pi = 3.14159265358979323846;
static constexpr float InvPi = 0.31830988618379067154;
static constexpr float Inv2Pi = 0.15915494309189533577;
static constexpr float Inv4Pi = 0.07957747154594766788;
static constexpr float PiOver2 = 1.57079632679489661923;
static constexpr float PiOver4 = 0.78539816339744830961;
static constexpr float Sqrt2 = 1.41421356237309504880;

// Global inline functions
template <typename T>
CIVET_CPU_GPU inline void swap(T &a, T &b) {
	T tmp = std::move(a);
	a = std::move(b);
	b = std::move(tmp);
}

CIVET_CPU_GPU
inline float lerp(float t, float v1, float v2) {
	return (1 - t) * v1 + t * v2;
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

} // namespace civet

#endif // CIVET_H