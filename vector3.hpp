
#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include "precision.hpp"

class Vector3{
private:

    real data[3];
    
public:

    static Vector3 ORIGIN;    

    Vector3(){
        data[0] = data[1] = data[2] = 0.0;
    }


    Vector3(real x, real y, real z){
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }


    Vector3(real x){
        data[0] = data[1] = data[2] = x;
    }


    real x() const{
        return data[0];
    }


    real y() const{
        return data[1];
    }


    real z() const{
        return data[2];
    }


    void setX(real x){
        data[0] = x;
    }


    void setY(real y){
        data[1] = y;
    }


    void setZ(real z){
        data[2] = z;
    }


    real lengthSquared() const {
        return data[0] * data[0] 
            + data[1] * data[1] 
            + data[2] * data[2];
    }


    real length() const {
        return sqrtReal(lengthSquared());
    }

    
    Vector3 operator-() const{
        Vector3 result(-data[0], -data[1], -data[2]);
        return result;
    }


    Vector3 operator+(const Vector3& v) const{
        Vector3 result(data[0] + v.data[0], 
            data[1] + v.data[1], 
            data[2] + v.data[2]);
        return result;
    }


    Vector3 operator-(const Vector3& v) const{
        Vector3 result(data[0] - v.data[0], 
            data[1] - v.data[1], 
            data[2] - v.data[2]);
        return result;
    }


    // Dot product
    real operator*(const Vector3& v) const{
        return data[0] * v.data[0]
            + data[1] * v.data[1]
            + data[2] * v.data[2];
    }


    // cross product
    Vector3 operator%(const Vector3& v) const{
        Vector3 result(
            data[1] * v.data[2] - data[2] * v.data[1],
            data[2] * v.data[0] - data[0] * v.data[2],
            data[0] * v.data[1] - data[1] * v.data[0] 
        );
        return result;
    }


    // Scalar product
    Vector3 operator*(real s) const{
        return Vector3(data[0] * s, data[1] * s, data[2] * s);
    }


    void normalize(){
        real d = length();
        if(std::abs(d) > EPSILON){
            data[0] /= d;
            data[1] /= d;
            data[2] /= d;
        }
    }


    Vector3 normalized()const{
        Vector3 result(*this);
        result.normalize();
        return result;
    }


    // Access using brackets
    const real& operator[](int index) const {
        assert((index >= 0 && index <= 2) && "Index out of bounds");
        return data[index];
    }

    
    real& operator[](int index) {
        assert((index >= 0 && index <= 2) && "Index out of bounds");
        return data[index];
    }


    // Checks if a function is in a certain cuboid bounding box
    bool inBounds(const std::pair<Vector3, Vector3>& bounds) const {
        for (int i = 0; i < 3; ++i) {
            if ((*this)[i] < bounds.first[i] || (*this)[i] > bounds.second[i]){
                return false;
            }
        }
        return true;
    }
    
};

#endif