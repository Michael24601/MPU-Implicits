// For convenience we will define a some vector functions that work 
// for nD vectors.

#ifndef VECTOR_HELPER_H
#define VECTOR_HELPER_H

#include "precision.hpp"
#include <stdexcept>
#include <vector>


real operator*(const std::vector<real>& u, 
    const std::vector<real>& v){

    assert((v.size() == u.size() && v.size() != 0)
        && "Vector dimensions must match for the dot product\n");  

    real result = 0;
    for(int i = 0; i < u.size(); i++){
        result += u[i] * v[i];
    }
    return result;
}


std::vector<real> operator+(const std::vector<real>& u, 
    const std::vector<real>& v){
        
    assert((v.size() == u.size() && v.size() != 0)
        && "Vector dimensions must match for addition\n");

    std::vector<real> result;
    result.resize(v.size());
    for(int i = 0; i < u.size(); i++){
        result[i] = u[i] + v[i];
    }
    return result;
}


std::vector<real> operator-(const std::vector<real>& u, 
    const std::vector<real>& v){
        
    assert((v.size() == u.size() && v.size() != 0)
        && "Vector dimensions must matchfor subtraction\n");  

    std::vector<real> result;
    result.resize(v.size());
    for(int i = 0; i < u.size(); i++){
        result[i] = u[i] - v[i];
    }
    return result;
}


std::vector<real> operator*(const std::vector<real>& u, 
    real scalar){

    std::vector<real> result;
    result.resize(u.size());
    for(int i = 0; i < u.size(); i++){
        result[i] = u[i] * scalar;
    }
    return result;
}


#endif
