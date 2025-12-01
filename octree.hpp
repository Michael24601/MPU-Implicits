
#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include "precision.hpp"
#include "vector3.hpp"
#include "kdTree3.hpp"


class Octree{
private:

    class Node{
    public:

        // Points associated with this node (storing only the index).
        std::vector<int> pointIndex;
        int depth;
        // Centroid of the cube associated with the node
        Vector3 centroid;
        real radius;

        Node(int depth, const Vector3& centroid) : 
            depth{depth}, centroid{centroid}{}

        
        real diagonal(){
            return halfPower[depth];
        }
    };


    // This only includes leafs, there is no advantage (as far
    // as this project goes) to storing the full tree.
    std::vector<Node> nodes;

public:

    Octree() {}

    
    // Given a cube centered at (0, 0, 0), with a main diagonal of size
    // 1, subdivides the cube such that some criteria is met.
    void subdivideSpace(int depth, const Vector3& center,
        const KdTree3& tree){

        real diagonal = halfPower[depth];
        real halfDiagonal= halfPower[depth+1];

        // We can assume that if we are at depth > 1,
        // that the parent of this node contained at least N_MIN points,
        // otherwise we would not have recursed.
        // So it is sufficient to search the paren't nodes for points 
        // contained inside this node.
        real radius = diagonal * ALPHA;
        std::vector<int> points;
        tree.rangeQuery(center, radius, points);

        // The empty ball case, no further subdivision
        if(points.size() == 0){
            
            return;
        }
        
        //
        while(points.size() < N_MIN){
            radius += LAMBDA;
            tree.rangeQuery(center, radius, points);
        }
        
        
    }


};

#endif