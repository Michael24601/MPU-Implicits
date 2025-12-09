
#ifndef BULK_INSERT_TEST_H
#define BULK_INSERT_TEST_H

#include "insert.hpp"
#include <random>


void bulkInsertTest(){

    std::vector<Point> points = generatePoints(1000, -10.0, 10.0, 987654);
    KdTree3 tree;

    tree.bulkInsert(points);
    
    std::vector<Point> subSet;
    tree.rangeQuery(Vector3::ORIGIN, 5.0, subSet);

    for(const Point& p: subSet){
        std::cout << (p.getPoint() - Vector3::ORIGIN).length() << "\n";
    }
}

#endif