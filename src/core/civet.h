#ifndef CIVET_H
#define CIVET_H

#include <iostream>
#include "gpu/testCuda.h"

#define CIVET_CPU_GPU __host__ __device__
#define CIVET_GPU_CODE __device__

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
template <typename T>
class Bounds2;
template <typename T>
class Bounds3;
class Ray;

class Medium;


// Global constants
#ifdef _MSC_VER
#define MaxFloat std::numeric_limits<float>::max()
#define Infinity std::numeric_limits<float>::infinity()
#else
static constexpr float MaxFloat = std::numeric_limits<float>::max();
static constexpr float Infinity = std::numeric_limits<float>::infinity();
#endif

// Global inline functions
inline float lerp(float t, float v1, float v2) { return (1 - t) * v1 + t * v2; }
}

#endif // CIVET_H