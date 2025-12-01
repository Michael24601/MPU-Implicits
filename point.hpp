
#ifndef POINT_H
#define POINT_H

#include "vector3.hpp"

// Oriented Point
class Point{
private:

    Vector3 point;
    Vector3 normal;

public:

    Point(real x, real y, real z, real nx, real ny, real nz)
        : point(x, y, z), normal(nx, ny, nz){
            real normalLength = normal.length();
            if(normalLength > 1.0 && abs(normalLength) >= EPSILON){
                normal = normal * (1.0 / normalLength);
            }
        }

};

#endif