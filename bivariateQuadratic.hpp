
#ifndef BIVARIATE_QUADRATIC_H
#define BIVARIATE_QUADRATIC_H

#include "localFitFunction.hpp"
#include "matrix.hpp"
#include "monomialBasis.hpp"
#include "conjugateGradient.hpp"
#include "weightFunction.hpp"
#include "vector2.hpp"

// Bivariate Quadratic polynomial
class BivariateQuadratic: public LocalFitFunction{
private:

    // The coefficients A, B, C, D, E, F    
    real coefficients[6];

    // The coordinate system of the bivariate is (u, v, w),
    // such that the plane (u, v) is orthogonal to the normal n,
    // and centered at the origin c (both c and n are expressed
    // in the original (x, y, z) coordinate system).
    // Also w is the distance along the normal.
    
    Vector3 normal;
    Vector3 origin;
    // basis
    Vector3 uBasis;
    Vector3 vBasis;


    // We call this before we fit the function, and it sets the basis
    // including the orthonormal 2D basis u, v, the normal, which is
    // also our 3rd basis vector w, and the origin.
    void setBasis(const Vector3& origin, const Vector3& normal){

        this->origin = origin;
        // (we normalize just in case, though it should already be normal)
        this->normal = normal.normalized();

        // First we generate u, an arbitrary vector on the plane
        // orthogonal to the normal.
        // We pick some random vector r not parallel to the normal
        // and get u from their cross product.
        Vector3 r(0, 0, 1);
        if(std::abs(std::abs(r * this->normal) - 1.0) < EPSILON){
            r = Vector3(1, 0, 0);
        }
        this->uBasis = (normal % r).normalized();
        
        // We then generate v such that it is perpendicular to u
        // and the normal.
        this->vBasis = (normal % uBasis).normalized();
    }


    // Transform a point from xyz to uvw.
    Vector3 transform(const Vector3& point) const {

        // This ensures the point is at the new origin
        Vector3 xyz = point - origin;
        Vector3 uvw;

        // We then project the point xyz, onto to the plane.
        uvw.setX(uBasis * xyz);
        uvw.setY(vBasis * xyz);

        // We then calculate the distance of xyz along the normal,
        // which is by definition the dot product.
        uvw.setZ(xyz * normal);

        return uvw;
    }


    // Returns the matrix M we get when we differentiate the
    // energy. See the report for the details on how I derived it.
    // As input we take the points p.
    // Also returns the vector b in Nc = d.
    // We do these at the same time to save time, as they share
    // certain computations.
    void getLinearSystem(
        const std::vector<Point>& points,
        const Vector3& center, real radius,
        Matrix& m, std::vector<real>& b
    ) const {
        
        // First we will compute the vector containing the
        // monomial basis function values for each point and auxiliary
        // point.
        // We then create the outer product matrix, and multiply it
        // by the weights.
        // And then we assemble M and b.

        // There are 6 basis functions in 2D with degree 2,
        // so the matrix is 6 by 6.
        Matrix sumP(6, 0.0);
        std::vector<real> vec(6, 0);
        for(int i = 0; i < points.size(); i++){

            // First we transform the points so that they are in uvw.
            Vector3 uvw = transform(points[i].getPoint());

            // There are 6 basis functions in 2D with degree 2.
            // (we ignore the w coordinate).
            Vector2 point(uvw.x(), uvw.y());
            std::vector<real> basis = MonomialBasis::evaluate(point, 2);
            // Outerproduct
            Matrix outerProduct(basis);   
                     
            // The b-spline (notice that we use the xyz coordinates
            // for this).
            real bSplineValue = WeightFunction::evaluate(center, 
                radius, points[i].getPoint());
            sumP = sumP + (outerProduct * bSplineValue);
            vec = Matrix::add(vec, Matrix::scale(basis, bSplineValue * uvw.z()));
        }

        // Finally we set the results
        m = sumP;
        b = vec;
    }


public:


    BivariateQuadratic(const Vector3& origin, const Vector3& normal): 
        coefficients{0,0,0,0,0,0}{

        setBasis(origin, normal);
    }


    // Fits the quadratic given a set of points.
    void fit(
        const std::vector<Point>& points,
        const Vector3& center, real radius
    ){
        Matrix m(6);
        std::vector<real> b(6);
        getLinearSystem(points, center, radius, m, b);
        // We just use a linear system solver
        std::vector<real> coefficients = ConjugateGradient::solve(m, b);
        // Now the basis functions are ordered such that
        // the coefficients correspond to A, C, B, D, E, F.
        this->coefficients[0] = coefficients[0];
        this->coefficients[1] = coefficients[2];
        this->coefficients[2] = coefficients[1];
        this->coefficients[3] = coefficients[3];
        this->coefficients[4] = coefficients[4];
        this->coefficients[5] = coefficients[5];
    }


    // The quadratic has a formula w - (Au^2 + 2Buv + Cv^2 + Du + Ev + F).a
    // Note however that the input is is given in the original (x, y, z)
    // coordinate system.
    real evaluate(const Vector3& input) const override{
        Vector3 uvw = transform(input);
        real u{uvw.x()}, v{uvw.y()}, w{uvw.z()};
        real result = w - (coefficients[0] * u * u +
            coefficients[1] * u * v * 2 + coefficients[2] * v * v +
            coefficients[3] * u + coefficients[4] * v + coefficients[5]);
        return result;
    }

    
    Vector3 evaluateGradient(const Vector3& input) const override{

        // To evaluate the gradient at a point xyz, we first
        // convert it to uvw, get the gradient nabla uvw, and
        // then we use the multivariate chain rule to get back to xyz.
        // In particular, we know that a linear transformation A
        // is a function f: R^3 -> R^3, so it will have a Jacobian,
        // which will happen to be the matrix A itself.
        // By the multivariate chain rule:
        // f(g(x, y, z))    = J_g^T * (nabla_{uvw} f(u, v, w))
        //                  = A^T * nabla_{uvw} f(u, v, w)
        // So we just need to multiply the gradient by the transpose
        // of the basis matrix we used to transform xyz to uvw.
        Vector3 uvw = transform(input);
        real u{uvw.x()}, v{uvw.y()}, w{uvw.z()};
        // The gradient (wrt to UVW):
        Vector3 gradient(
            - (2 * coefficients[0] * u + 2 * coefficients[1] * v 
                + coefficients[3]),
            - (2 * coefficients[1] * u + 2 * coefficients[2] * v 
                + coefficients[4]), 
            1.0
        );

        // Multiplying by the transpose means the basis vectors,
        // which were rows, become columns. So we take a
        // linear combination instead of using the dot product.
        Vector3 gradientXYZ(uBasis * gradient.x()
            + vBasis * gradient.y() + normal * gradient.z());
        return gradientXYZ;
    }

};

#endif