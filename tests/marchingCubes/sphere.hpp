
#ifndef SPHERE_TEST_H
#define SPHERE_TEST_H

#include "../../marchingCubes.hpp"

// The sphere is centered at (0.5, 0.5, 0.5) with radius 0.25.
real sphereSignedDistance(const Vector3& p){
    // The -radius is here for the sign
    return (p - Vector3(0.5)).length() - 0.25;
}

// Tests marching cube using a sphere implicit function
void sphereTest(){

    real (*f)(const Vector3&) = sphereSignedDistance;
    // Uses 20 cubes per direction
    marchingCubes(f, Vector3(1.0), 0, 30, 30, 30, false,
        "output/sphere_test.obj");
}

#endif