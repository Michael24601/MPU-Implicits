
#ifndef TIME_COMPLEXITY_TEST_H
#define TIME_COMPLEXITY_TEST_H

#include "insert.hpp"
#include <random>
#include <chrono>

void timeComplexityTest(){

    for(int i = 10'000; i <= 1'000'000; i += 10'000){
        std::vector<Point> points = generatePoints(i, -100.0, 100.0, 987321);
        KdTree3 tree;
        tree.bulkInsert(points);

        auto start = std::chrono::high_resolution_clock::now();

        // Subset with roughly half the points
        std::vector<Point> subSet;
        tree.rangeQuery(Vector3::ORIGIN, 5.0, subSet);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << elapsed.count() << "\n";
    }
}

#endif