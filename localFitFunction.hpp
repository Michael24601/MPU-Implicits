
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
};


#endif