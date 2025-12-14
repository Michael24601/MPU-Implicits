
#ifndef POINT_H
#define POINT_H

#include "vector3.hpp"

// Oriented Point
class Point{
private:

    Vector3 point;
    Vector3 normal;

public:

    Point() : point(0), normal(0){}

    Point(const Vector3& pos, const Vector3& n) 
        : point(pos), normal(n) {}

    Point(real x, real y, real z, real nx, real ny, real nz)
        : point(x, y, z), normal(nx, ny, nz){
        real normalLength = normal.length();
        if(normalLength > 1.0 && std::abs(normalLength) >= EPSILON){
            normal = normal * (1.0 / normalLength);
        }
    }


    const Vector3& getPoint() const{
        return point;
    }


    const Vector3& getNormal() const{
        return normal;
    }


    void setPoint(const Vector3& point){
        this->point = point;
    }

};

#endif