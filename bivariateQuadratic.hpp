
#ifndef BIVARIATE_QUADRATIC_H
#define BIVARIATE_QUADRATIC_H

#include "localFitFunction.hpp"


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
        this->normal = normal;

        // First we generate u, an arbitrary vector on the plane
        // orthogonal to the normal.
        // We pick some random vector r not parallel to the normal
        // and get u from their cross product.
        Vector3 r(0, 0, 1);
        if(abs(r * normal - 1.0) < EPSILON){
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


public:


    BivariateQuadratic(const Vector3& origin, const Vector3& normal): 
        coefficients{0,0,0,0,0,0}{

        setBasis(origin, normal);
    }


    // Fits the quadric given a set of points (also initializes
    // the normal and center)
    void fit(const std::vector<Point>& point){
        // TODO:
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
        // then multiply each element by the basis vectorx u, v, and w 
        // to get the gradient of xyz (because the gradient is a normal,
        // that is, it is covariant, so transforming the gradient
        // from uvw to xyz uses the same matrix as transforming
        // a contravariant vector from xyz to uvw, that is).
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
        Vector3 gradientXYZ(
            uBasis * gradient,
            vBasis * gradient,
            normal * gradient
        );
        return gradientXYZ;
    }

};

#endif