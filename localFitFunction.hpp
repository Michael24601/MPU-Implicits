
#ifndef LOCAL_FIT_FUNCTION_H
#define LOCAL_FIT_FUNCTION_H

#include <vector>
#include "precision.hpp"
#include "vector3.hpp"
#include "point.hpp"


// Pure virtual class for a local fir function Q(x)
class LocalFitFunction{
public:
    virtual real evaluate(const Vector3& input) const = 0;
    virtual Vector3 evaluateGradient(const Vector3& input) const = 0;

    // This calculates the max-norm approximation error,
    // used to decide if we subdivide further.
    // We assume the given points are the one inside the sphere
    // of the octree to whom this local function belongs.
    real approximationError(const std::vector<Point>& points) const{
        real max = 0;
        for(int i = 0; i < points.size(); i++){
            // Note, we use std::abs, since abs apparently causes
            // floating point numbers to become 0.
            real error = std::abs(evaluate(points[i].getPoint()));
            error /= evaluateGradient(points[i].getPoint()).length();
            if(error > max){
                max = error;
            }
        }
        return max;
    }
};


#endif