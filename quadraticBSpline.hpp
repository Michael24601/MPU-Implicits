
#ifndef QUADRATIC_B_SPLINE_H
#define QUADRATIC_B_SPLINE_H

#include "vector3.hpp"

class QuadraticBSpline{

private:

    // Evalutes the quadratic B-spline on the local support
    // [0, 3]. See the report for the full derivation.
    static inline constexpr real evaluateBSpline(real t){
        if(t < 0.0 || t > 3.0){
            return 0;
        }
        else if(t <= 1.0){
            return t * t * 0.5;
        }
        else if(t <= 2.0){
            return (-2.0 * t * t + 6.0 * t - 3.0) * 0.5;
        }
        else{
            return (3.0 - t) * (3.0 - t) * 0.5;
        }
    }

public:

    static real evaluate(const Vector3& center, real radius, 
        const Vector3& p){
        
        // The input is ||x - c||/R
        real input = (p - center).length() / radius;
        
        // First we map the input, which is naturally in [0, 1]
        // to be in [0, 3/2], then we can evaluate the B-spline.
        return evaluateBSpline(input * 1.5);
    }

};

#endif