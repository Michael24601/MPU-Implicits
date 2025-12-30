
#ifndef OBJECT_TEST_H
#define OBJECT_TEST_H

#include "../../octree.hpp"
#include "../../marchingCubes.hpp"

#include <cmath>
#include <algorithm>
#include <random>


struct Triangle {
    Point p0, p1, p2;
    float area;
};

float rand01() {
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX); 
}

// Helper: compute area of a triangle
float triangleArea(const Vector3& a, const Vector3& b, const Vector3& c) {
    return ((b - a) % (c - a)).length() * 0.5f; // cross product
}

// Sample a point uniformly on a triangle using barycentric coordinates
Point sampleTriangle(const Triangle& tri) {
    float u = rand01(); // random float [0,1)
    float v = rand01();
    if (u + v > 1.0f) { u = 1.0f - u; v = 1.0f - v; }
    float w = 1.0f - u - v;

    Vector3 pos = tri.p0.getPoint() * w 
        + tri.p1.getPoint() * u + tri.p2.getPoint() * v;
    Vector3 n   = (tri.p0.getNormal() * w + tri.p1.getNormal() 
        * u + tri.p2.getNormal() * v).normalized();
    return {pos, n};
}

std::vector<Point> sampleObjTriangles(const std::string& filename, int k) {
    std::ifstream file(filename);
    if (!file.is_open()) return {};

    std::vector<Vector3> vertices;
    std::vector<Vector3> vertexNormals;
    std::vector<Triangle> triangles;
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {
            float x, y, z; ss >> x >> y >> z;
            vertices.push_back(Vector3(x, y, z));
        } else if (prefix == "vn") {
            float x, y, z; ss >> x >> y >> z;
            vertexNormals.push_back(Vector3(x, y, z));
        } else if (prefix == "f") {
            int v[3], vn[3] = {-1, -1, -1};
            for (int i = 0; i < 3; i++) {
                std::string vertStr; ss >> vertStr;
                size_t firstSlash = vertStr.find('/');
                size_t lastSlash = vertStr.rfind('/');
                if (firstSlash != std::string::npos) {
                    v[i] = std::stoi(vertStr.substr(0, firstSlash)) - 1;
                    if (lastSlash != firstSlash)
                        vn[i] = std::stoi(vertStr.substr(lastSlash + 1)) - 1;
                } else {
                    v[i] = std::stoi(vertStr) - 1;
                }
            }

            Vector3& v0 = vertices[v[0]];
            Vector3& v1 = vertices[v[1]];
            Vector3& v2 = vertices[v[2]];

            Vector3 fallback = ((v1 - v0) % (v2 - v0)).normalized();
            Vector3 n0 = (vn[0] >= 0) ? vertexNormals[vn[0]] : fallback;
            Vector3 n1 = (vn[1] >= 0) ? vertexNormals[vn[1]] : fallback;
            Vector3 n2 = (vn[2] >= 0) ? vertexNormals[vn[2]] : fallback;

            Triangle tri;
            tri.p0 = {v0, n0};
            tri.p1 = {v1, n1};
            tri.p2 = {v2, n2};
            tri.area = triangleArea(v0, v1, v2);

            triangles.push_back(tri);
        }
    }

    // Build PDF from triangle areas
    std::vector<float> cdf(triangles.size());
    float totalArea = 0;
    for (size_t i = 0; i < triangles.size(); i++) {
        totalArea += triangles[i].area;
        cdf[i] = totalArea;
    }

    std::vector<Point> points;
    for (int i = 0; i < k; i++) {
        float r = ((float)rand() / RAND_MAX) * totalArea;

        // binary search on CDF
        auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
        int triIdx = std::distance(cdf.begin(), it);

        points.push_back(sampleTriangle(triangles[triIdx]));
    }

    return points;
}


void scalePoints(std::vector<Point>& points) {
    
    // Bounding box
    Vector3 minPt = points[0].getPoint();
    Vector3 maxPt = points[0].getPoint();
    for (const auto& p : points) {
        minPt.setX(std::min(minPt.x(), p.getPoint().x()));
        minPt.setY(std::min(minPt.y(), p.getPoint().y()));
        minPt.setZ(std::min(minPt.z(), p.getPoint().z()));
        maxPt.setX(std::max(maxPt.x(), p.getPoint().x()));
        maxPt.setY(std::max(maxPt.y(), p.getPoint().y()));
        maxPt.setZ(std::max(maxPt.z(), p.getPoint().z()));
    }

    // The largest extent
    float maxExtent = std::max({maxPt.x() - minPt.x(), maxPt.y() - 
        minPt.y(), maxPt.z() - minPt.z()});
    // Full range must be smaller than (-0.5/sqrt(3), 0.5/sqrt(3))
    float targetExtent = INV_SQRT3 * 0.9;
    float scale = targetExtent / maxExtent;

    for (auto& p : points) {
        p.setPoint((p.getPoint() - (minPt + maxPt) * 0.5f) * scale);
    }
}


// This function just writes to an obj. It only writes oriented points,
// no faces.
void writeToFile(const std::vector<Point>& points, 
    const std::string& filename) {

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    for (const Point& p : points) {
        file << "v " << p.getPoint().x() << " " 
            << p.getPoint().y()  << " " << p.getPoint().z()  << "\n";
    }
    for (const Point& p : points) {
        file << "vn " << p.getNormal().x() << " " 
            << p.getNormal().y()  << " " << p.getNormal().z()  << "\n";
    }
    for (int i = 0; i < points.size(); i++) {
        int k = i+1;
        file << "f " << k << "//" << k
            << " " << k << "//" << k
            << " " << k << "//" << k << "\n";
    }

    file.close();
}


void objectTest(){

    // We first generate the points
    std::vector<Point> points = sampleObjTriangles("input/bunny.obj", 15000);
    scalePoints(points);
    KdTree3 tree;
    tree.bulkInsert(points);

    // Then we create the octree
    Octree octree;
    octree.subdivideOctree(tree);

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
        0, 100, 100, 100, "output/bunny_test.obj");

    writeToFile(points, "output/bunny_point_cloud.obj");
}

#endif