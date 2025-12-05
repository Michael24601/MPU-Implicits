
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


    // Transform a point from xyz to uvw.
    Vector3 transform(const Vector3& point) const {

        // This ensures the resulting point is centered at the
        // given origin.
        Vector3 xyz = point - origin;
        Vector3 uvw;

        // First we generate u, an arbitrary vector on the plane
        // orthogonal to the normal.
        // We pick some random vector r not parallel to the normal
        // and get u from their cross product.
        Vector3 r(0, 0, 1);
        if(abs(r * normal - 1.0) < EPSILON){
            r = Vector3(1, 0, 0);
        }
        Vector3 u = (normal % r).normalized();
        
        // We then generate v such that it is perpendicular to u
        // and the normal.
        Vector3 v = (normal % u).normalized();

        // We then project the point xyz, onto to the plane.
        uvw.setX(u * xyz);
        uvw.setY(v * xyz);

        // We then calculate the distance of xyz along the normal,
        // which is by definition the dot product.
        uvw.setZ(xyz * normal);

        return uvw;
    }


public:


    BivariateQuadratic(){}


    // Fits the quadric given a set of points (also initializes
    // the normal and center)
    void fit(const std::vector<Point>& point, const Vector3& normal,
        const Vector3& center){
        
        
    }


    // The quadratic has a formula w - (Au^2 + 2Buv + Cv^2 + Du + Ev + F).
    // Note however that the input is is given in the original (x, y, z)
    // coordinate system.
    real evaluate(const Vector3& input) const override{
        Vector3 uvw = transform(input);
        real u{uvw.x()}, v{uvw.y()}, w{uvw.z()};
        real result = w - (coefficients[0] * u * u +
            coefficients[1] * u * v * 2 +
            coefficients[2] * v * v +
            coefficients[3] * u +
            coefficients[4] * v +
            coefficients[5]);
        return result;
    }

};

#endif