
#ifndef RAY_MARCHING_SPHERE_TEST_H
#define RAY_MARCHING_SPHERE_TEST_H

#include "../../rayMarching.hpp"

// The sphere is centered at (0.5, 0.5, 0.5) with radius 0.25.
real sphereSignedDistanceRay(const Vector3& p){
    // The -radius is here for the sign
    return (p - Vector3(0.5)).length() - 0.25;
}

Vector3 color(const Vector3& p){
    return Vector3(255.0);
}

// Tests marching cube using a sphere implicit function
void raySphereTest(){

    real (*f)(const Vector3&) = sphereSignedDistanceRay;
    Vector3 (*colorF)(const Vector3&) = color;

    RayMarching::run(
        f, colorF, Vector3(0.25, 0, 0),
        Vector3(0, 0.25, 0), Vector3(0.5, 0.5, -1.0), 1.0,
        Vector3(0.5), SQRT3, 512, 512, Vector3(0.0), 
        "output/ray_sphere_test.png"
    );
}

#endif