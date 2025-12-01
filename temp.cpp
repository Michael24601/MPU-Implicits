
#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include "precision.hpp"
#include "vector3.hpp"

// The minimum number of points in a leaf node
constexpr int N_MIN = 15;

class Octree{
private:

    class Node{
    public:

        // Points associated with this node (storing only the index).
        // Only applies to leaves.
        std::vector<int> pointIndex;
        int index;
        bool isLeaf;

        Node(int index) : index{index}, isLeaf{true}{}
    };

    // The current trees max depth (root is at 0)
    int depth;
    std::vector<Node> nodes;

public:

    Octree() : depth{0} {
        // Root
        nodes.push_back(Node(0));
    }


    // Size of an octree with some depth (root is at depth 0).
    inline int size(int depth){
        return (std::pow(8, depth+1) - 1) / 7;
    }


    int incrementDepth(){
        depth++;
        nodes.resize(size(depth));
    }


    Node& getNode(int parentIndex, int childIndex){
        int index = parentIndex * 8 + childIndex;
        if(index >= 0 && index < nodes.size()){
            return nodes[index];
        }
        else{
            throw std::invalid_argument("Node index is invalid\n");
        }
    }


    int getNodeDepth(int index) const{
        return static_cast<int>(
            std::floor(std::log(7 * index + 1) / std::log(8))
        );
    }

    
    // Subdivides the tree at this particular node
    void subdivideNode(int index){

        // If the current node is at the max possible depth
        // (this is faster than explicitely calculating the node depth).
        if(index * 8 + 1 >= nodes.size()){
            incrementDepth();
        }

        nodes[index].isLeaf = false;
        
    }


};

#endif