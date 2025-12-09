
#ifndef GENERAL_QUADRIC_H
#define GENERAL_QUADRIC_H

#include "localFitFunction.hpp"


// General Quadric
class GeneralQuadric: public LocalFitFunction{
private:

    // Elements of the matrix A (lower left triangle)
    real matrix[6];

    // The b vector
    Vector3 vec;

    // The c scalar offset
    real scalar;

public:


    GeneralQuadric() : matrix{0,0,0,0,0,0}, vec(0,0,0), scalar(0) {}


    // Fits the quadric given a set of points
    void fit(const std::vector<Point>& point, Vector3* corners){
        // TODO:
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


    // Since the matrix A is symmetric, the gradient is 2Ax + b
    Vector3 evaluateGradient(const Vector3& input) const override{
        Vector3 result(
            2 * (matrix[0] * input.x() + matrix[1] * input.y() 
                + matrix[3] * input.z()) + vec.x(),
            2 * (matrix[1] * input.x() + matrix[2] * input.y() 
                + matrix[4] * input.z()) + vec.y(),
            2 * (matrix[3] * input.x() + matrix[4] * input.y() 
                + matrix[5] * input.z()) + vec.z()
        );
        return result;
    }

};


#endif