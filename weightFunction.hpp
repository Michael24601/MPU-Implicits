
#ifndef WEIGHT_FUNCTION_H
#define WEIGHT_FUNCTION_H

#include "bSpline.hpp"
#include "inverseDistanceWeight.hpp"

class WeightFunction{

    public:

    // Helper function, so that we can easily swicth weight functions
    // without needing to change every call to the weight function.
    static real evaluate(const Vector3& center, real radius, const Vector3&p){
        return QuadraticBSpline::evaluate(center, radius, p);
    }
};

#endif