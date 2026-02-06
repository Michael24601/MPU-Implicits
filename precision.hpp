// File that every other file in the project should include.
// Just defines the precision used (float or double)
// as well as some useful constants.

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

// Stored here for faster lookup (1/2^n)
// Must ensure it is at least MAX_DEPTH+2 in length
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