
#ifndef RAY_MARCHING_OBJECT_TEST_H
#define RAY_MARCHING_OBJECT_TEST_H

#include "../../rayMarching.hpp"
#include "../evaluation/object.hpp"


Vector3 objectColor(const Vector3& p){
    return Vector3(155.0);
}

// Tests marching cube using a sphere implicit function
void rayObjectTest(){

    // We first generate the points
    std::vector<Point> points = sampleObjTriangles("input/bunny.obj", 1);
    KdTree3 tree;
    tree.bulkInsert(points);

    // Then we create the octree
    Octree octree;
    octree.subdivideOctree(tree);

    // Function Q
    auto lambdaFunction = [&octree](const Vector3& p) -> real {
        return octree.evaluate(p);
    };

    auto lambdaColor = [&octree](const Vector3& p) -> Vector3 {
        int depth = octree.getDepth(p);
        // Blue is shallow
        real blue = 255.0 * (MAX_DEPTH - depth) / MAX_DEPTH;
        // Red is deep
        real red = 255.0 * depth / MAX_DEPTH;
        return Vector3(red, 0, blue);
    };

    Vector3 (*colorF)(const Vector3&) = objectColor;

    RayMarching::run(
        lambdaFunction, lambdaColor, Vector3(0.25, 0, 0),
        Vector3(0, 0.25, 0), Vector3(0, 0, 1.0), -1.0,
        Vector3(0.5), SQRT3, 1024, 1024, Vector3(0.0), 
        "output/ray_bunny_test.png"
    );
}

#endif