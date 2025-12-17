
#ifndef PRECISION_H
#define PRECISION_H

#include <cmath>
#include <iostream>
#include <cassert>

// Precision
typedef double real;
inline real (*sqrtReal)(real) = ::sqrt;

// Constants
constexpr real PI = 3.141592653589793;
constexpr real HALF_PI = 1.570796326794897;
constexpr real INV_SQRT3 = 0.5773502691896257;
constexpr real SQRT3 = 1.7320508075688772;

// Epsilon used to check if floating point numbers are 0
constexpr real EPSILON = 1e-8;

// The minimum number of points in a leaf node
constexpr int N_MIN = 15;

// The factor used in the radius computation
constexpr real ALPHA = 0.75;

// The step size when incrementing the radius
constexpr real LAMBDA = 0.1;

// The maximum depth the octree can go
constexpr int MAX_DEPTH = 8;

// The maximum accepted approximation error
// (Must be careful with this one, too small a value and the octree
// never stops subsdividing and the max depth is always reached)
constexpr real EPSILON_ZERO = 0.01;

// The conjugate gradient parameters
constexpr real MAX_CG_ERROR = 1e-8;
constexpr real MAX_CG_ITERATIONS = 100;

// Controls ray marcher step size (will be scaled by the scene)
// This needs to be small to avoid staircase artifcats, where
// large portions of the rendered image look like a solid color.
constexpr real RAY_MARCHER_STEP = 0.001;
constexpr real FINITE_DIFFERENCE_STEP = 0.001;

// Stored here for faster lookup (1/2^n)
// Ensure it is at least MAX_DEPTH+2 in length
constexpr real halfPower[22] = {
    1.0,
    0.5,
    0.25,
    0.125,
    0.0625,
    0.03125,
    0.015625,
    0.0078125,
    0.00390625,
    0.001953125,
    0.0009765625,
    0.00048828125,
    0.000244140625,
    0.0001220703125,
    0.00006103515625,
    0.000030517578125,
    0.0000152587890625,
    0.00000762939453125,
    0.000003814697265625,
    0.0000019073486328125,
    0.00000095367431640625,
    0.000000476837158203125};

#endif