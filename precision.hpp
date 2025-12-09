
#ifndef PRECISION_H
#define PRECISION_H

#include <cmath>

// Precision
typedef float real;
inline float (*sqrtReal)(float) = ::sqrtf;

// Constants
constexpr real PI = 3.141592653589793;
constexpr real HALF_PI = 1.570796326794897;
constexpr real INV_SQRT3 = 0.5773502691896257;

// Epsilon used to check if floating point numbers are 0
constexpr real EPSILON = 1e-8;

// The minimum number of points in a leaf node
constexpr int N_MIN = 15;

// The factor used in the radius computation
constexpr real ALPHA = 0.75;

// The step size when incrementing the radius
constexpr real LAMBDA = 0.1;

// The maximum depth the octree can go
constexpr int MAX_DEPTH = 20;

// The maximum accepted approximation error
constexpr real EPSILON_ZERO = 0.5;

// The conjugate gradient parameters
constexpr real MAX_CG_ERROR = 1e-06;
constexpr real MAX_CG_ITERATIONS = 50;

// Stored here for faster lookup (1/2^n)
constexpr real halfPower[MAX_DEPTH + 2] = {
    1.0f,
    0.5f,
    0.25f,
    0.125f,
    0.0625f,
    0.03125f,
    0.015625f,
    0.0078125f,
    0.00390625f,
    0.001953125f,
    0.0009765625f,
    0.00048828125f,
    0.000244140625f,
    0.0001220703125f,
    0.00006103515625f,
    0.000030517578125f,
    0.0000152587890625f,
    0.00000762939453125f,
    0.000003814697265625f,
    0.0000019073486328125f,
    0.00000095367431640625f,
    0.000000476837158203125f
};

#endif