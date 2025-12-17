
#ifndef INVERSE_DISTANCE_WEIGHT_H
#define INVERSE_DISTANCE_WEIGHT_H

#include "vector3.hpp"

class InverseDistanceWeight{

private:

    // Function that returns a if a > 0, 0 otherwise
    static inline constexpr real positivePart(real a){
        return std::max(a, 0.0);
    }

    // This is just the inverse distance function istelf
    static inline constexpr real evaluateInverseDistance(real radius, 
        real distance){
    
        real value = positivePart(radius - distance) / (radius * distance);
        return value * value;
    }   

public:

    static real evaluate(const Vector3& center, real radius, 
        const Vector3& p){
        
        real distance = (center - p).length();
        return evaluateInverseDistance(radius, distance);
    }

};

#endif