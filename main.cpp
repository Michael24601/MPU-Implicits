
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <chrono>

#include "octree.hpp"
#include "marchingCubes.hpp"
#include "rayMarching.hpp"
#include "meshUtil.hpp"


// Parameters

real Octree::ALPHA = 0.75;
real Octree::LAMBDA = 0.1;
real Octree::EPSILON_ZERO = 0.001;
unsigned int Octree::N_MIN = 15;
unsigned int Octree::MAX_DEPTH = 8;
bool Octree::USE_PIECEWISE_POLYNOMIALS = false;

real RayMarching::RAY_MARCHER_STEP = 0.001;
real RayMarching::FINITE_DIFFERENCE_STEP = 0.001;

real ConjugateGradient::MAX_CG_ERROR = 1e-8;
unsigned int ConjugateGradient::MAX_CG_ITERATIONS = 30;

bool MarchingCubes::USE_MIDDLE_POINT = false;


// Which output we want (possibly all)
constexpr bool OUTPUT_MARCHING_CUBES = true;
constexpr bool OUTPUT_RAY_MARCHING = false;
constexpr bool OUTPUT_POINT_CLOUD = false;

int main(){

    // We first generate the points, by sampling randomly from
    // a mesh. We choose a triangle at random (using a pdf
    // weighted by the surface area of each triangle), and then
    // sample randomly on the triangle a point using barycentric
    // coordinates.
    std::vector<Point> points = sampleObjTriangles("input/bunny.obj", 30000);

    // We then scale the points such that they are in a cube centered
    // at the origin with a main diagonal of length 1.
    scalePoints(points);

    // Then we insert the points into a kdTree.
    KdTree3 tree;
    tree.bulkInsert(points);

    auto start = std::chrono::high_resolution_clock::now();

    // Then we create the octree, and let the fitting procedure complete
    Octree octree;
    octree.subdivideOctree(tree);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Time elapsed: " << elapsed.count() << "\n";

    // We can also do marching cubes, using this lambda
    // as an input density.
    auto lambdaFunction = [&octree](const Vector3& p) -> real {
        // Note, we flip the sign as marching cubes expects a higher
        // density inside.
        // The 0.5 * INV_SQRT3 offset is just because
        // marching cubes expects the shape between (0, 0, 0)
        // and some other point "range".
        return octree.evaluate((p - Vector3(0.5 * INV_SQRT3))) * -1.0;
    };


    if(OUTPUT_MARCHING_CUBES){
        MarchingCubes::run(lambdaFunction, Vector3(INV_SQRT3), 
            0, 200, 200, 200, "output/bunny_depth_test.obj");
    }


    // We can also use ray marching

    // We define a color function that shows depth at each point
    auto lambdaDepthColor = [&octree](const Vector3& p) -> Vector3 {

        int depth = octree.getDepth(p);
        // Maps depth into is in [0, 1]
        real t = real(depth) / real(Octree::MAX_DEPTH);

        // This just makes the scaling non-linear, since most
        // points on the surface are in deeper levels.
        real power = 2; 
        real tn = std::pow(t, power);

        // Red is deep, blue is shallow, and we interpolate
        // between the two using the scaled depth
        real red  = 255.0 * tn;
        real blue = 255.0 * (1.0 - tn);
        return Vector3(red, 0, blue);
    };

    auto lambdaGreyColor = [](const Vector3& p) -> Vector3 {
        return Vector3(155.0);
    };

    auto rayMarchingLambdaFunction = [&octree](const Vector3& p) -> real {
        return octree.evaluate(p);
    };


    if(OUTPUT_RAY_MARCHING){
        RayMarching::run(
            rayMarchingLambdaFunction, lambdaGreyColor, 
            Vector3(0.25, 0, 0), Vector3(0, 0.25, 0), 
            Vector3(0, 0, 1.0), -1.0, 
            Vector3(0.5), SQRT3, 1024, 1024, Vector3(0.0), 
            "output/ray_bunny_test.png"
        );
    }
    
    // We can also output the point cloud for reference
    if(OUTPUT_POINT_CLOUD){
        writeToFile(points, "output/bunny_point_cloud.obj");
    }

    return 0;
}