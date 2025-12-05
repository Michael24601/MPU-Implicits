
#ifndef KD_TREE_H
#define KD_TREE_H

#include "vector3.hpp"
#include "point.hpp"
#include <vector>

// A 3D kd-tree implementation, used to store 3D points
class KdTree3{
private:

    class Node{
    public:

        Node* left;
        Node* right;
        Point point;


        Node(const Point& point) : point{point}, 
            left{nullptr}, right{nullptr}{}

    };

    typedef Node* NodePtr;
    NodePtr root;

    void destructorAux(NodePtr ptr){
        if(ptr){
            if(ptr->left) destructorAux(ptr->left);
            if(ptr->right) destructorAux(ptr->right);
            delete ptr;
        }
    }

    void rangeQueryAux(const Vector3& center, real radius, 
        std::vector<Point>& points, NodePtr ptr, int axis) const {

        if(!ptr){
            return;
        }

        Vector3 difference = center - ptr->point.getPoint();

        if(difference.lengthSquared() <= radius * radius){
            points.push_back(ptr->point);
        }

        // If the point is within the range with respect to this
        // axis only, we recurse on that side.
        if(difference[axis] <= radius){
            rangeQueryAux(center, radius, points, 
                ptr->left, (axis+1)%3);
        }
        if(difference[axis] >= -radius){
            rangeQueryAux(center, radius, points, 
                ptr->right, (axis+1)%3);
        }
    }

public:

    KdTree3() : root{nullptr}{}

    ~KdTree3(){
        destructorAux(root);
    }


    void insert(const Point& point){

        if(!root){
            root = new Node(point);
            return;
        }

        NodePtr ptr = root;
        int axis = 0;

        // It's like a binary tree but at each node we alternate 
        // which axis we compare.
        while(ptr){
            if(ptr->point.getPoint()[axis] > point.getPoint()[axis]){
                if(ptr->left){
                    ptr = ptr->left;
                }
                else{
                    ptr->left = new Node(point);
                    return;
                }
            }
            else{
                if(ptr->right){
                    ptr = ptr->right;
                }
                else{
                    ptr->right = new Node(point);
                    return;
                }
            }
            axis = (axis + 1) % 3;
        }
    }

    
    // Returns all points within a certain range (returns only the indexes)
    void rangeQuery(const Vector3& center, real radius, 
        std::vector<Point>& points) const {

        rangeQueryAux(center, radius, points, root, 0);
    }

};


#endif