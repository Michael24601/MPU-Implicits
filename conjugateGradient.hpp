// A conjugate gradient solver.
// Note that I used https://en.wikipedia.org/wiki/Conjugate_gradient_method
// as a reference.

#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "matrix.hpp"
#include "vectorHelper.hpp"

class ConjugateGradient{

public:

    // This assumes that A is an SPD matrix
    static std::vector<real> solve(const Matrix& A, 
        const std::vector<real>& b){

        assert(A.rows() == b.size() 
            && "Matrix and vector dimensions must match\n");

        typedef std::vector<real> Vector;

        const real MAX_ERROR_SQUARED = MAX_CG_ERROR * MAX_CG_ERROR;

        Vector x0(A.columns(), 0);
        Vector r0 = b - (A * x0);
        Vector p0 = r0;

        if(r0 * r0 < MAX_ERROR_SQUARED){
            return x0;
        }

        int iteration = 0;
        // Though we assume A is SPD, we have a max iteration guard
        // in case of no convergence
        while(iteration < MAX_CG_ITERATIONS){

            real dotProduct = r0 * r0;
            real alpha = dotProduct;
            Vector prod = A * p0;
            alpha /= p0 * prod;
            Vector x1 = x0 + (p0 * alpha);
            Vector r1 = r0 - (prod * alpha);

            real beta = r1 * r1;
            if(beta < MAX_ERROR_SQUARED){
                x0 = x1;
                break;
            }

            beta /= dotProduct;
            p0 = r1 + (p0 * beta);
            x0 = x1;
            r0 = r1;
            iteration++;
        }

        return x0;
    }
};

#endif