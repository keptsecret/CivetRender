#ifndef CIVET_H
#define CIVET_H

#include <iostream>
#include "gpu/testCuda.h"

#define CIVET_CPU_GPU __host__ __device__
#define CIVET_GPU_CODE __device__

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

#endif // CIVET_H