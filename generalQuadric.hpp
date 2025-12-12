
#ifndef GENERAL_QUADRIC_H
#define GENERAL_QUADRIC_H

#include "localFitFunction.hpp"
#include "matrix.hpp"
#include "monomialBasis.hpp"
#include "quadraticBSpline.hpp"
#include "conjugateGradient.hpp"
#include "kdTree3.hpp"

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
    static void getLinearSystem(
        const std::vector<Point>& points, 
        const std::vector<Vector3>& aux, const std::vector<real>& diag,
        const Vector3& center, real radius,
        Matrix& m, std::vector<real>& b
    ){
        
        // First we will compute the vector containing the
        // monomial basis function values for each point and auxiliary
        // point.
        // We then create the outer product matrix, and multiply it
        // by the weights.
        // And then we assemble M and b.

        // There are 10 basis functions in 3D with degree 2,
        // so the matrix is 10 by 10.
        Matrix sumP(10, 0.0);
        real sumBSplines = 0;
        for(int i = 0; i < points.size(); i++){
            // There are 10 basis functions in 3D with degree 2.
            std::vector<real> basis = 
                MonomialBasis::evaluate(points[i].getPoint(), 2);
            // Outerproduct
            Matrix outerProduct(basis);            
            // The b-spline
            real bSplineValue = QuadraticBSpline::evaluate(center, 
                radius, points[i].getPoint());
            sumP = sumP + (outerProduct * bSplineValue);
            sumBSplines += bSplineValue;
        }
        // We then divide the matrix sums by the sum of b-splines
        sumP = sumP * (1.0 / sumBSplines);
        
        // Then we do the same thing with the auxiliary points,
        // using them to generate the second matrix term and the
        // vector b.
        Matrix sumQ(10, 0.0);
        std::vector<real> vec(10, 0.0);
        for(int i = 0; i < aux.size(); i++){
            std::vector<real> basis = MonomialBasis::evaluate(aux[i], 2);
            // Outerproduct
            Matrix outerProduct(basis);
            sumQ = sumQ + outerProduct;
            vec = Matrix::add(vec, Matrix::scale(basis, diag[i])); 
        }
        sumQ = sumQ * (1.0 / aux.size());

        // Finally we set the results
        m = sumP + sumQ;
        b = vec;
    }


    // Filters which aux points are usable and which ones are not.
    // Returns both the points and their diagonals
    static void filterAuxiliaryPoints(
        const std::vector<Point>& points,
        const std::vector<Vector3>& aux, int k,
        std::vector<Vector3>& usableAux,
        std::vector<real>& diags
    ){
        // For each aux point q, we will find the k nearest neighbors,
        // and compute n_i . (q - p_i) with its k neighbors.
        // It is usable if the sign matches.
        KdTree3 tree;
        // We have to copy the points since the bulk insert modifies them
        std::vector<Point> copy(points);
        tree.bulkInsert(copy);

        for(int i = 0; i < aux.size(); i++){

            std::vector<Point> closest;
            tree.kNNSearch(aux[i], k, closest);
            real prod;
            bool flagSame = true;
            // For the diagonals
            real sum = 0.0;

            for(int j = 0; j < closest.size(); j++){
                real dotProduct = closest[j].getNormal() 
                    * (aux[i] - closest[j].getPoint());
                sum += dotProduct;
                if(j == 0){
                    prod = dotProduct;
                }
                else if(dotProduct * prod < 0){
                    flagSame = false;
                    break;
                }
            }

            if(flagSame) {
                usableAux.push_back(aux[i]);
                diags.push_back(sum / closest.size());
            }
        }
    }


public:


    GeneralQuadric() : matrix{0,0,0,0,0,0}, vec(0,0,0), scalar(0) {}


    // Fits the quadric given a set of points and auxiliary points.
    // Returns false if for any reason the fitting doesn't work
    // and the cell needs to be subdivided
    bool fit(const std::vector<Point>& points, 
        const std::vector<Vector3>& aux,
        const Vector3& center, real radius
    ){

        if(points.empty() || aux.empty()){
            throw std::invalid_argument("Can't fit quadric on empty points\n");
        }

        // First we check which auxiliary points are actually usable
        std::vector<Vector3> usableAux;
        std::vector<real> diags;
        filterAuxiliaryPoints(points, aux, 6, usableAux, diags);

        // If no aux points are usable, we can't fit the function
        if(usableAux.empty()){
            return false;
        }

        // The matrix is 10 by 10
        Matrix m(10);
        std::vector<real> b;
        getLinearSystem(points, usableAux, diags, center, radius, m, b);
        // We just use a linear system solver
        std::vector<real> coefficients = ConjugateGradient::solve(m, b);

        matrix[0] = coefficients[0];
        matrix[2] = coefficients[1];
        matrix[5] = coefficients[2];
        matrix[1] = coefficients[3];
        matrix[3] = coefficients[4];
        matrix[4] = coefficients[5];

        vec.setX(coefficients[6]);
        vec.setY(coefficients[7]);
        vec.setZ(coefficients[8]);
        scalar = coefficients[9];

        return true;
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