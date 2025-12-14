
#ifndef EVALUATION_H
#define EVALUATION_H

#include "../../octree.hpp"
#include "../../kdTree3.hpp"
#include "../../marchingCubes.hpp"


std::vector<Point> generatePointsOnSphere(int n, int m){
    std::vector<Point> points(n*m);
    for(int i = 0; i < n; i++){
        real theta = PI * (i + 0.5) / n;
        for(int j = 0; j < m; j++){
            real phi = 2.0 * PI * j / m;

            real x = std::sin(theta) * std::cos(phi);
            real y = std::sin(theta) * std::sin(phi);
            real z = std::cos(theta);

            Vector3 pos(x, y, z);
            // We want the main diagonal of the cube to be 1, so we need
            // the sphere to fit inside that. We do that by multiplying
            // the pos by 1/2 * sqrt3.
            // We then multiply by 0.9 to ensure it fits well.
            pos = pos * (0.5 * INV_SQRT3 * 0.9);
            Vector3 normal = pos.normalized();
            points[i*m + j] = Point(pos, normal);
        }
    }
    
    return points;
}


real sphereSignedDistanceReference(const Vector3& p){
    // The -radius is here for the sign
    return (p - Vector3(0.5 * INV_SQRT3)).length() - 0.5 * INV_SQRT3 * 0.9;
}


void evaluationTest(){

    // We first generate the points
    std::vector<Point> points = generatePointsOnSphere(100, 100);
    KdTree3 tree;
    tree.bulkInsert(points);

    // Then we create the octree
    Octree octree;
    octree.subdivideOctree(tree);

    // Then we evaluate it at different places to check if the result
    // is as expected.
    // We check that the surface points are near 0.
    // We can't check the signed distance function as there is no
    // guarantee that we get an approximation of it, just that we have
    // non-zero values that are positive outside, and negative inside.
    int n = 10, m = 10;
    real error = 0;
    for(int i = 0; i < n; i++){
        real theta = PI * (i + 0.5) / n;
        for(int j = 0; j < m; j++){
            real phi = 2.0 * PI * j / m;

            real x = std::sin(theta) * std::cos(phi);
            real y = std::sin(theta) * std::sin(phi);
            real z = std::cos(theta);

            Vector3 pos(x, y, z);
            pos = pos * 0.5 * INV_SQRT3 * 0.9;
            error += octree.evaluate(pos);
        }
    }

    error /= (n*m);
    std::cout << error << "\n";

    // We can also do marching cubes
    auto lambdaFunction = [&octree](const Vector3& p) -> real {
        // Note, we flip the sign as marching cubes expects a higher
        // density inside.
        // The 0.5 * INV_SQRT3 offset is just because
        // marching cubes expects the shape between (0, 0, 0)
        // and some other point "range".
        return octree.evaluate((p - Vector3(0.5 * INV_SQRT3))) * -1.0;
    };

    MarchingCubes::run(lambdaFunction, Vector3(INV_SQRT3), 
        0, 50, 50, 50, false, "output/sphere_evaluation_test.obj");

    // We can run it on the actual signed distance function for
    // comparison
    MarchingCubes::run(sphereSignedDistanceReference, Vector3(INV_SQRT3), 
        0, 50, 50, 50, false, "output/sphere_evaluation_reference.obj");

}


#endif