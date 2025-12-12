// While writing this class, I mostly cross-referenced my code
// with the wikipedia article on kd-Trees:
// https://en.wikipedia.org/wiki/K-d_tree

#ifndef KD_TREE_H
#define KD_TREE_H

#include <vector>
#include <algorithm>
#include <queue>
#include "vector3.hpp"
#include "point.hpp"

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


    // Helper class
    struct Neighbor{
        Point point;
        real distanceSquared;

        // Required for the priority queue
        bool operator<(const Neighbor& n) const {
            return distanceSquared < n.distanceSquared;
        }
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


    // Finds median point in given range (both sides inclusive), 
    // inserts it into the tree, and then recurses on either side of it.
    // The axis alternates per step.
    NodePtr bulkInsertAux(std::vector<Point>& points, int start, int end, 
        int axis){

        // The node at which this function was called must be
        // a leaf if start > end.
        if(start > end){
            return nullptr;
        }
        
        // We find the middle index
        int middle = (start+end)/2;

        // We then use this function, which sorts only the middle element
        // of the array, essentially ensuring all elements to the
        // left are smaller, and all elements to the right are larger
        // (within the given range only, with respect to the current axis).
        std::nth_element(
            points.begin() + start,
            points.begin() + middle,
            points.begin() + end + 1,
            [axis](const Point& p1, const Point& p2) { 
                return p1.getPoint()[axis] < p2.getPoint()[axis]; 
            }
        );

        // We then insert the middle element into the tree as a new node.
        NodePtr node = new Node(points[middle]);
        // We then recurse
        node->left = bulkInsertAux(points, start, middle-1, (axis+1) % 3);
        node->right = bulkInsertAux(points, middle+1, end, (axis+1) % 3);
        
        return node;
    }


    // This is a nearest neighbor search expedited by the kd tree.
    void kNNSearchAux(const Vector3& p, int k, 
        NodePtr ptr, int axis, std::priority_queue<Neighbor>& knn) const {
        
        if(!ptr){
            return;
        }

        // First we check if the current point is closer to the reference
        // point than the furthest point in the priority queue.
        // Unless the queue has less than k points, in which case we just
        // add it.
        real lengthSquared = (p - ptr->point.getPoint()).lengthSquared();
        if(knn.size() < k){
            knn.push(Neighbor{ptr->point, lengthSquared});
        }
        else if(lengthSquared < knn.top().distanceSquared){
            knn.pop();
            knn.push(Neighbor{ptr->point, lengthSquared});
        } 

        // Next we check which subtree is closer so we can recurse
        // into it.
        // We check the distance between p and the current point
        // with respect to the current axis.
        real distance = p[axis] - ptr->point.getPoint()[axis];
        NodePtr closer, farther;
        if(distance < 0){
            closer = ptr->left;
            farther = ptr->right;
        }
        else{
            closer = ptr->right;
            farther = ptr->left;
        }

        kNNSearchAux(p, k, closer, (axis+1)%3, knn);

        // We only search the farther subtree if we need to, that is,
        // if the distance to the point (with respect to the current
        // axis) is smaller than the largest distance to a point
        // in the priority queue
        if(distance * distance < knn.top().distanceSquared){
            kNNSearchAux(p, k, farther, (axis+1)%3, knn);
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

    
    // Returns all points within a certain range.
    void rangeQuery(const Vector3& center, real radius, 
        std::vector<Point>& points) const {

        rangeQueryAux(center, radius, points, root, 0);
    }


    // Inserts bulk data all at once into the tree.
    // Balanced kd-trees perform better, and since we have access
    // to all the data at once, instead of a self-balancing tree,
    // we can ensure balance by inserting the data in a specific way.
    // We insert the median, then recursively search for the median
    // in the two halves that we have left, alternating the axis
    // we use to find the median.
    // This also removes the need to use the insert function.
    // Note that this assumes the tree is emtpy.
    // Also note that the points vector will be modified (some of its
    // elements will be sorted).
    void bulkInsert(std::vector<Point>& points){
        if(root){
            return;
        }
        else{
            root = bulkInsertAux(points, 0, points.size()-1, 0);
        }
    }


    // This is a k-nearest neighbor search.
    // It can be used to for the quadric fitting.
    void kNNSearch(const Vector3& p, int k, 
        std::vector<Point>& points) const {

        std::priority_queue<Neighbor> knn;

        kNNSearchAux(p, k, root, 0, knn);

        // Then we just add the result into the vector
        points.reserve(points.size() + knn.size());

        while(!knn.empty()) {
            points.push_back(knn.top().point);
            knn.pop();
        }
    }

};


#endif