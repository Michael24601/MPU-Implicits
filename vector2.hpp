
#ifndef VECTOR2_H
#define VECTOR2_H

#include <iostream>
#include "precision.hpp"

class Vector2{
private:

    real data[2];
    
public:

    static Vector2 ORIGIN;    

    Vector2(){
        data[0] = data[1] = 0.0;
    }


    Vector2(real x, real y){
        data[0] = x;
        data[1] = y;
    }


    Vector2(real x){
        data[0] = data[1] = x;
    }


    real x() const{
        return data[0];
    }


    real y() const{
        return data[1];
    }


    // Access using brackets
    const real& operator[](int index) const {
        if(index < 0 || index > 1){
            throw std::invalid_argument("Index out of bounds");
        }
        else{
            return data[index];
        }
    }
};

#endif