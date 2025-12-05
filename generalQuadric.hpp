
#ifndef GENERAL_QUADRIC_H
#define GENERAL_QUADRIC_H

#include "localFitFunction.hpp"


// General Quadric
class GeneralQuadric: public LocalFitFunction{
private:

    // Elements of the matrix A (lowe right triangle)
    real matrix[6];

    // The b vector
    Vector3 vec;

    // The c scalar offset
    real scalar;

public:


    GeneralQuadric(){}


    // Fits the quadric given a set of points
    void fit(const std::vector<Point>& point, Vector3* corners){
        
    }


    // The quadric has a formula x^TAx + b^Tx + c
    real evaluate(const Vector3& input) const override{
        real result = scalar;
        result += vec * input;
        result += matrix[0] * input.x() * input.x();
        result += 2 * matrix[1] * input.x() * input.y();
        result += matrix[2] * input.y() * input.y();
        result += 2 * matrix[3] * input.x() * input.z();
        result += 2 * matrix[4] * input.y() * input.z();
        result += matrix[5] * input.z() * input.z(); 
        return result;
    }

};


#endif