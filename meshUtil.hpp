// File that contains some helper functions for reading and writing
// and sampling meshes.

#ifndef MESH_UTIL_H
#define MESH_UTIL_H


#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include "octree.hpp"


struct Triangle {
    Point p0, p1, p2;
    float area;
};


// Returns a random number between 0 and 1
float randomNumber() {
    static std::random_device device;
    static std::mt19937 r(device());
    static std::uniform_real_distribution<real> dis(0.0, 1.0);
    return dis(r);
}


// Returns area of triangle, using the cross product
float triangleArea(const Vector3& a, const Vector3& b, const Vector3& c) {
    return ((b - a) % (c - a)).length() * 0.5f;
}


// Samples a point uniformly on a triangle using barycentric coordinates
Point sampleTriangle(const Triangle& t) {

    real u = randomNumber();
    real v = randomNumber();
    // If outside of the triangle, we fold them back in
    if (u + v > 1.0) { 
        u = 1.0f - u; 
        v = 1.0 - v; 
    }
    float w = 1.0f - u - v;

    Vector3 position = t.p0.getPoint() * w 
        + t.p1.getPoint() * u + t.p2.getPoint() * v;
    Vector3 normal = (t.p0.getNormal() * w + t.p1.getNormal() * u 
        + t.p2.getNormal() * v).normalized();

    return Point(position, normal);
}


// Chooses triangles at random and samples them
std::vector<Point> sampleObjTriangles(const std::string& filename, int k) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file\n";
        return std::vector<Point>{};
    }


    // First we store the triangles from the mesh file
    std::vector<Vector3> vertices;
    std::vector<Vector3> vertexNormals;
    std::vector<Triangle> triangles;
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string prefix;
        ss >> prefix;

        // Vertex case
        if (prefix == "v") {
            real x, y, z; 
            ss >> x >> y >> z;
            vertices.push_back(Vector3(x, y, z));
        } 
        // Normal case
        else if (prefix == "vn") {
            real x, y, z; 
            ss >> x >> y >> z;
            vertexNormals.push_back(Vector3(x, y, z));
        } 
        // Face case
        else if (prefix == "f") {

            // We always expect to find vertices
            int v[3];
            // We use -1 to indicate that no normal was found
            int vn[3] = {-1, -1, -1};

            for (int i = 0; i < 3; i++) {
                std::string vertexString; 
                ss >> vertexString;
                int firstSlash = vertexString.find('/');
                int lastSlash = vertexString.rfind('/');

                // We know that the face format is
                // f v//n v//n v//n
                // Or
                // f v v v
                // Where v and n are the indexes of faces
                // and normals, starting at 1

                // If there is a first slash
                if (firstSlash != std::string::npos) {
                    // Then we record the vertex (number before slash)
                    // and subtract 1 since the arrays are 0 not 1 indexed.
                    v[i] = std::stoi(vertexString.substr(0, firstSlash)) - 1;

                    // If there is a second slash, we also have a normal
                    // (If we had one slash, the second number would
                    // not be a normal)
                    if (lastSlash != firstSlash){
                        vn[i] = std::stoi(vertexString.substr(lastSlash + 1)) 
                            - 1;
                    }
                }
                // Otherwise if we have no first slash, and the only
                // number we have is the vertex index.
                else {
                    v[i] = std::stoi(vertexString) - 1;
                }
            }

            Vector3 v0 = vertices[v[0]];
            Vector3 v1 = vertices[v[1]];
            Vector3 v2 = vertices[v[2]];

            // If for some reason the mesh does not contain any
            // vertex normals, we can use the geometric normal of
            // the triangle.
            Vector3 geoNormal = ((v1 - v0) % (v2 - v0)).normalized();
            Vector3 n0 = (vn[0] >= 0) ? vertexNormals[vn[0]] : geoNormal;
            Vector3 n1 = (vn[1] >= 0) ? vertexNormals[vn[1]] : geoNormal;
            Vector3 n2 = (vn[2] >= 0) ? vertexNormals[vn[2]] : geoNormal;

            Triangle triangle;
            triangle.p0 = Point(v0, n0);
            triangle.p1 = Point(v1, n1);
            triangle.p2 = Point(v2, n2);
            triangle.area = triangleArea(v0, v1, v2);

            triangles.push_back(triangle);
        }
    }

    // Once we have all triangles, we build a PDF using their
    // areas.


    // This contains the cumulative surface area,
    // that is, pdf[0] = area of triangle 0
    // pdf[1] = area of triangle 0 and 1 ...
    std::vector<real> pdf(triangles.size());
    float totalArea = 0;
    for (int i = 0; i < triangles.size(); i++) {
        totalArea += triangles[i].area;
        pdf[i] = totalArea;
    }

    // Then, to sample points, what we do is select randomly
    // a number between 0 and 1, and map it between 0 and total area.
    // Then wherever it ends up in the pdf array based
    // on the interval pdf[i-1] and pdf[i] it is in, tells us
    // which triangle we choose. This ensures the probability
    // of landing on one triangle is proportional to its surface
    // area relative to the total area, since ending up in
    // the "section" between pdf[i-1] and pdf[i] is exactly the
    // area.
    std::vector<Point> points;
    for (int i = 0; i < k; i++) {
        real random = randomNumber() * totalArea;
        auto it = std::lower_bound(pdf.begin(), pdf.end(), random);
        int index = std::distance(pdf.begin(), it);

        // We then sample the triangle at the specific index,
        // also randomly.
        points.push_back(sampleTriangle(triangles[index]));
    }

    return points;
}


// Scales points
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
    real maxExtent = std::max(maxPt.x() - minPt.x(), 
        std::max(maxPt.y() - minPt.y(), maxPt.z() - minPt.z()));
    // Full range must be smaller than (-0.5/sqrt(3), 0.5/sqrt(3)),
    // and we use the largest extent to determine the scaling factor
    // to get to that target.
    real targetExtent = INV_SQRT3 * 0.9;
    real scale = targetExtent / maxExtent;

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
        // The index for the vertices and normals starts at 1,
        // so we add 1.
        int k = i+1;
        file << "f " << k << "//" << k
            << " " << k << "//" << k
            << " " << k << "//" << k << "\n";
    }

    file.close();
}

#endif