
#ifndef MONOMIAL_BASIS_H
#define MONOMIAL_BASIS_H

#include <vector>
#include "vector3.hpp"
#include "vector2.hpp"

class MonomialBasis{

public:

    // Returns the monomial basis function evaluations in 3D
    // for a degree n polynomial.
    // The order in which we return them tells us the order in which
    // the coefficients are laid out in the c vector in the linear
    // system we have to solve for the optimization problem.
    // We only need up to degree 2 however.
    static std::vector<real> evaluate(const Vector3& p, int degree){
        switch(degree){
        case 0: return std::vector<real>{1};
        case 1: return std::vector<real>{p.x(), p.y(), p.z(), 1};
        case 2: return std::vector<real>{p.x()*p.x(), p.y()*p.y(), 
            p.z()*p.z(), 2*p.x()*p.y(), 2*p.x()*p.z(), 2*p.y()*p.z(),
            p.x(), p.y(), p.z(), 1};
        default: return std::vector<real>(0);
        }
    }


     // Returns the monomial basis function evaluations in 2D
    // for a degree n polynomial.
    static std::vector<real> evaluate(const Vector2& p, int degree){
        switch(degree){
        case 0: return std::vector<real>{1};
        case 1: return std::vector<real>{p.x(), p.y(), 1};
        case 2: return std::vector<real>{p.x()*p.x(), p.y()*p.y(), 
            2*p.x()*p.y(), p.x(), p.y(), 1};
        default: return std::vector<real>(0);
        }
    }


};

#endif