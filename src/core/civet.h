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
#include <iostream>
#include <iterator>
#include <new>
#include <string>
#include <vector>
#ifdef CIVET_HAVE_MALLOC_H
#include <malloc.h> // for _alloca, memalign
#endif
#ifdef CIVET_HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <thread>
#include <type_traits>
#include <typeinfo>
#include <utility>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <imgui.h>
#define IMGUI_DEFINE_MATH_OPERATORS
#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl3.h>
#include <imgui/imgui_internal.h>
#include <imgui/misc/cpp/imgui_stdlib.h>

#define ALLOCA(TYPE, COUNT) (TYPE*)alloca((COUNT) * sizeof(TYPE))

namespace civet {

// Global forward declarations
class Scene;
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
class Matrix4;

class Ray;
class RayDifferential;
struct Interaction;
class SurfaceInteraction;

class Shape;
class Primitive;
class GeometricPrimitive;
class TransformedPrimitve;
class GLModel;
class Skybox;

class Material;
template <typename T>
class Texture;
class BxDF;
class BSDF;
class BSSRDF;

class Medium;
class MediumInteraction;
struct MediumInterface;

class Light;
class VisibilityTester;
class AreaLight;
class GLLight;
class GLDirectionalLight;
class GLPointLight;

template <int n_spectrum_samples>
class CoefficientSpectrum;
class RGBSpectrum;
class SampledSpectrum;
///< Choose between spectrum types by uncommenting the relevant line
typedef RGBSpectrum Spectrum;
// typedef SampledSpectrum Spectrum;

class Camera;
struct CameraSample;
class ProjectiveCamera;
class Film;
class FilmTile;
class Filter;
class Sampler;

class MemoryArena;
template <typename T, int logBlockSize = 2>
class BlockedArray;
class RNG;

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
inline uint32_t floatToBits(float f) {
	uint32_t ui;
	memcpy(&ui, &f, sizeof(float));
	return ui;
}

inline float bitsToFloat(uint32_t ui) {
	float f;
	memcpy(&f, &ui, sizeof(uint32_t));
	return f;
}

inline uint64_t floatToBits(double f) {
	uint64_t ui;
	memcpy(&ui, &f, sizeof(double));
	return ui;
}

inline double bitsToFloat(uint64_t ui) {
	double f;
	memcpy(&f, &ui, sizeof(uint64_t));
	return f;
}

template <typename T>
CIVET_CPU_GPU inline void swapElem(T& a, T& b) {
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
		swapElem(t0, t1);
	}
	return true;
}

CIVET_CPU_GPU
inline float gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline float gammaCorrect(float value) {
	if (value <= 0.0031308f) {
		return 12.92f * value;
	}
	return 1.055f * std::pow(value, (float)(1.f / 2.4f)) - 0.055f;
}

inline float inverseGammaCorrect(float value) {
	if (value <= 0.04045f) {
		return value * 1.f / 12.92f;
	}
	return std::pow((value + 0.055f) * 1.f / 1.055f, (float)2.4f);
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
	return (Pi / 180.0f) * deg;
}

CIVET_CPU_GPU
inline float degrees(float rad) {
	return (180.0f / Pi) * rad;
}

inline float log2(float x) {
	const double inv_log2 = 1.442695040888963387004650940071;
	return std::log(x) * inv_log2;
}

inline int log2Int(uint32_t v) {
#if defined(CIVET_IS_MSVC)
	unsigned long lz = 0;
	if (_BitScanReverse(&lz, v)) return lz;
	return 0;
#else
	return 31 - __builtin_clz(v);
#endif
}

inline int log2Int(int32_t v) { return log2Int((uint32_t)v); }

inline int log2Int(uint64_t v) {
#if defined(CIVET_IS_MSVC)
	unsigned long lz = 0;
#if defined(_WIN64)
	_BitScanReverse64(&lz, v);
#else
	if  (_BitScanReverse(&lz, v >> 32))
		lz += 32;
	else
		_BitScanReverse(&lz, v & 0xffffffff);
#endif // _WIN64
	return lz;
#else  // CIVET_IS_MSVC
	return 63 - __builtin_clzll(v);
#endif
}

inline int log2Int(int64_t v) { return log2Int((uint64_t)v); }

template <typename Predicate>
CIVET_CPU_GPU int findInterval(int size, const Predicate& pred) {
	int first = 0, len = size;
	while (len > 0) {
		int half = len >> 1, middle = first + half;
		if (pred(middle)) {
			first = middle + 1;
			len -= half + 1;
		} else {
			len = half;
		}
	}
	return clamp(first - 1, 0, size - 2);
}

template <typename T>
inline CIVET_CONSTEXPR bool isPowerOf2(T v) {
	return v && !(v & (v - 1));
}

inline int32_t roundUpPow2(int32_t v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

inline int64_t roundUpPow2(int64_t v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v |= v >> 32;
	return v + 1;
}

inline float erfInv(float x) {
	float w, p;
	x = clamp(x, -.99999f, .99999f);
	w = -std::log((1 - x) * (1 + x));
	if (w < 5) {
		w = w - 2.5f;
		p = 2.81022636e-08f;
		p = 3.43273939e-07f + p * w;
		p = -3.5233877e-06f + p * w;
		p = -4.39150654e-06f + p * w;
		p = 0.00021858087f + p * w;
		p = -0.00125372503f + p * w;
		p = -0.00417768164f + p * w;
		p = 0.246640727f + p * w;
		p = 1.50140941f + p * w;
	} else {
		w = std::sqrt(w) - 3;
		p = -0.000200214257f;
		p = 0.000100950558f + p * w;
		p = 0.00134934322f + p * w;
		p = -0.00367342844f + p * w;
		p = 0.00573950773f + p * w;
		p = -0.0076224613f + p * w;
		p = 0.00943887047f + p * w;
		p = 1.00167406f + p * w;
		p = 2.83297682f + p * w;
	}
	return p * x;
}

inline float erf(float x) {
	// constants
	float a1 = 0.254829592f;
	float a2 = -0.284496736f;
	float a3 = 1.421413741f;
	float a4 = -1.453152027f;
	float a5 = 1.061405429f;
	float p = 0.3275911f;

	// Save the sign of x
	int sign = 1;
	if (x < 0) {
		sign = -1;
	}
	x = std::abs(x);

	// A&S formula 7.1.26
	float t = 1 / (1 + p * x);
	float y =
			1 -
			(((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(-x * x);

	return sign * y;
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

inline GLenum glCheckError(const char* func_def) {
	GLenum error;
	while ((error = glGetError()) != GL_NO_ERROR) {
		std::string error_message;
		switch (error) {
			case GL_INVALID_ENUM:
				error_message = "INVALID_ENUM";
				break;
			case GL_INVALID_VALUE:
				error_message = "INVALID_VALUE";
				break;
			case GL_INVALID_OPERATION:
				error_message = "INVALID_OPERATION";
				break;
			case GL_OUT_OF_MEMORY:
				error_message = "OUT_OF_MEMORY";
				break;
			case GL_INVALID_FRAMEBUFFER_OPERATION:
				error_message = "INVALID_FRAMEBUFFER_OPERATION";
				break;
			default:
				error_message = "UNKNOWN_ERROR";
				break;
		}
		std::cerr << func_def << " " << error << " " << error_message << '\n';
	}
	return error;
}

} // namespace civet

#endif // CIVET_H