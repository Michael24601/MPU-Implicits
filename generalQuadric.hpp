
#ifndef GENERAL_QUADRIC_H
#define GENERAL_QUADRIC_H

#include "localFitFunction.hpp"
#include "matrix.hpp"
#include "monomialBasis.hpp"
#include "quadraticBSpline.hpp"

// General Quadric
class GeneralQuadric: public LocalFitFunction{
private:

    // Elements of the matrix A (lower left triangle)
    real matrix[6];

    // The b vector
    Vector3 vec;

    // The c scalar offset
    real scalar;


    // Returns the matrix M we get when we differentiate the
    // energy. See the report for the details on how I derived it.
    // As input we take the points p, the subset q,
    // and the diagonals d.
    // Also returns the vector b in Mc = b.
    // We do these at the same time to save time, as they share
    // certain computations.
    void getLinearSystem(
        const std::vector<Point>& points, 
        const std::vector<Vector3>& aux, const std::vector<real>& diag,
        const Vector3& center, real radius,
        Matrix& m, std::vector<real>& b
    ) const{
        
        // First we will compute the vector containing the
        // monomial basis function values for each point and auxiliary
        // point.
        // We then create the outer product matrix, and multiply it
        // by the 

        // For the points P

        // There are 10 basis functions in 3D with degree 2,
        // so the matrix is 10 by 10.
        Matrix sumP(10);
        real sumBSplines = 0;
        for(int i = 0; i < points.size(); i++){
            // There are 10 basis functions in 3D with degree 2.
            std::vector<real> basis = 
                MonomialBasis::evaluate(points[i].getPoint(), 2);
            // Outerproduct
            Matrix outerProduct = Matrix(basis);
            // The b-spline
            real bSplineValue = QuadraticBSpline::evaluate(center, 
                radius, points[i].getPoint());
            sumP = sumP + outerProduct * bSplineValue;
            sumBSplines += bSplineValue;
        }
        // We then divide the matrix sums by the sum of b-splines
        sumP = sumP * (1.0 / sumBSplines);
        
        // Then we do the same thing with the auxiliary points,
        // using them to generate the second matrix term and the
        // vector b.
        Matrix sumQ(10);
        std::vector<real> vec;
        for(int i = 0; i < aux.size(); i++){
            std::vector<real> basis = MonomialBasis::evaluate(aux[i], 2);
            // Outerproduct
            Matrix outerProduct = Matrix(basis);
            sumQ = sumQ + outerProduct;
            vec = Matrix::add(vec, Matrix::scale(basis, diag[i])); 
        }
        sumQ = sumQ * (1.0 / aux.size());

        // Finally we set the results
        m = sumP + sumQ;
        b = vec;
    }

public:


    GeneralQuadric() : matrix{0,0,0,0,0,0}, vec(0,0,0), scalar(0) {}


    // Fits the quadric given a set of points and auxiliary points
    void fit(const std::vector<Point>& point, 
        const std::vector<Vector3>& aux,
        const Vector3& center, real radius
    ){
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