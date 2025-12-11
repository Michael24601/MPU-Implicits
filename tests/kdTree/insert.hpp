
#ifndef INSERT_TEST_H
#define INSERT_TEST_H

#include "../../kdTree3.hpp"
#include <random>

std::vector<Point> generatePoints(int n, float min, float max, int seed){
    std::vector<Point> points;
    points.reserve(n);

    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> distribution(min, max);

    for (size_t i = 0; i < n; i++) {
        float x = distribution(rng);
        float y = distribution(rng);
        float z = distribution(rng);

        // normal doesn't matter for the purpose of this test
        points.push_back(Point(x, y, z, 0, 0, 1));
    }

    return points;
}

void insertTest(){

    std::vector<Point> points = generatePoints(1000, -10.0, 10.0, 123456);
    KdTree3 tree;

    for(const Point& point: points){
        tree.insert(point);
    }

    std::vector<Point> subSet;
    tree.rangeQuery(Vector3::ORIGIN, 5.0, subSet);

    for(const Point& p: subSet){
        std::cout << (p.getPoint() - Vector3::ORIGIN).length() << "\n";
    }
}

#endif