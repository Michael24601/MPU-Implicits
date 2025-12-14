
#ifndef BUNNY_TEST_H
#define BUNNY_TEST_H

#include "../../octree.hpp"
#include "../../marchingCubes.hpp"

#include <cmath>
#include <algorithm>
#include <random>


void randomBarycentric(real& u, real& v) {
    std::random_device rd;
    static std::mt19937 mt(rd());
    static std::uniform_real_distribution<real> distribution(0.0, 1.0);
    u = distribution(mt);
    v = distribution(mt);
    if (u + v > 1.0) {
        u = 1 - u;
        v = 1 - v;
    }
}


std::vector<Point> sampleObjTriangles(const std::string& filename, int k) {
    std::ifstream file(filename);
    if (!file.is_open()) return {};

    std::vector<Vector3> vertices;
    std::vector<Vector3> vertexNormals;
    std::vector<Point> points;
    std::string line;

    int triangleIndex = 0;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {
            float x, y, z;
            ss >> x >> y >> z;
            vertices.push_back(Vector3(x, y, z));
        } else if (prefix == "vn") {
            float x, y, z;
            ss >> x >> y >> z;
            vertexNormals.push_back(Vector3(x, y, z));
        } else if (prefix == "f") {
            int v[3], vn[3] = {-1, -1, -1};
            for (int i = 0; i < 3; i++) {
                std::string vertStr;
                ss >> vertStr;

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


            for (int i = 0; i < k; i++) {
                real u, vCoord;
                randomBarycentric(u, vCoord);
                real w = 1.0f - u - vCoord;

                Vector3 samplePoint = v0 * w + v1 * u + v2 * vCoord;
                Vector3 sampleNormal = (n0 * w + n1 * u + n2 * vCoord).normalized();

                points.push_back({samplePoint, sampleNormal});
            }

            triangleIndex++;
        }
    }

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

    return points;
}


void bunnyTest(){

    // We first generate the points
    std::vector<Point> points = sampleObjTriangles("input/bunny.obj", 10);
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
        0, 100, 100, 100, false, "output/bunny_test.obj");
}

#endif