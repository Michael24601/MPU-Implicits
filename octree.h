
#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include "precision.hpp"
#include "vector3.hpp"
#include "localFitFunction.hpp"
#include "kdTree3.hpp"


// Cube corners
const Vector3 cubeCorners[8] = {
    Vector3(-1, -1, -1),
    Vector3(1, -1, -1),
    Vector3(-1, 1, -1),
    Vector3(-1, -1, 1),
    Vector3(1, 1, -1),
    Vector3(1, -1, 1),
    Vector3(-1, 1, 1),
    Vector3(1, 1, 1)
};


class Octree{
private:

    // Note that the nodes store the localFitFunction, but not the
    // points contained in the sphere.
    class Node{
    private:

        int depth;

        // Centroid and radius of the cube associated with the node
        Vector3 centroid;
        real radius;
        // The local fit function
        LocalFitFunction* localFitFunction;
        
        // The 8 children
        Node* child[8] = {nullptr};
        bool leaf;

    public:

        Node(int depth, const Vector3& centroid) : 
            depth{depth}, centroid{centroid}, leaf{false},
            child{nullptr, nullptr, nullptr, nullptr, nullptr, 
                nullptr, nullptr, nullptr}{

            // By default the radius is ALPHA * diagonal()    
            radius = ALPHA * getDiagonal();
         }

        
        ~Node(){
            delete localFitFunction;
        }

        
        int getDepth() const{
            return depth;
        }
        

        Vector3 getCentroid() const{
            return centroid;
        }


        real getRadius() const{
            return radius;
        }


        void incrementRadius(){
            // By default we always add LAMBDA when incrementing
            // the radius.
            radius += LAMBDA;
        }


        bool isLeaf() const{
            return leaf;
        }


        void setRadius(real radius){
            this->radius = radius;
        }


        // The diagonal of each subdivision is some power (1/2)^n.
        real getDiagonal(){
            return halfPower[depth];
        }


        // Subdivides at this node by initializing 8 children
        void subdivide(){
            if(depth == MAX_DEPTH){
                throw std::exception("Tree is at maximum depth");
                return;
            }

            if(leaf){
                leaf = false;
                // The diagonal for this cube is halfPower[depth],
                // so the offset diagonal is at the quarter, that is is,
                // halfPower[depth+2], which we know exists since we
                // haven't reached the max depth.
                // We multiply iy by 1/sqrt3 to get the offset in x,
                // y and z.
                real offset = halfPower[depth+2] * INV_SQRT3;

                for(int i = 0; i < 8; i++){
                    child[i] = new Node(depth+1, 
                        centroid + (cubeCorners[i] * offset));
                }
            }
        }

    
        // Assuming that this node contains the given points,
        // this function fits a local function to the points.
        // Note that this function assumes that points is not empty.
        void fitLocalFunction(const std::vector<Point>& points){

            // It first calculates the maximum angle theta between
            // the normals of the points and the normalized arithmetic
            // mean of the points.
            Vector3 normalMean = std::accumulate(
                points.begin(), points.end(), 0.0, 
                [](Vector3& acc, const Point& s) {
                    return acc + s.getNormal();
                }
            );
            normalMean.normalize();

            auto it = std::max_element(
                points.begin(), points.end(), 
                [&normalMean](const Point& p1, const Point& p2) {
                    return p1.getNormal() * normalMean 
                        < p2.getNormal() * normalMean;
                }
            );
            real theta = (it == points.end() ? 0 : 
                it->getNormal() * normalMean);

            if(points.size() > 2 * N_MIN && theta > HALF_PI){

                // General quadric case
                GeneralQuadric* quadric = new GeneralQuadric();
                
                Vector3 corners[8];
                // We need half of this cube's diagonal
                real offset = halfPower[depth+1] * INV_SQRT3;

                for(int i = 0; i < 8; i++){
                    corners[i] = centroid + (cubeCorners[i] * offset);
                }

                quadric->fit(points, corners);
                this->localFitFunction = quadric;
            }
            else if(points.size() > 2 * N_MIN){
                
                // Bivariate
            }
            else{
                // Edges and corners
            }

        }

    };


    typedef Node* NodePtr;
    NodePtr root;


    // Given a node representing a cube in the octree, this function
    // computes the local fit function, and if some criteria is met,
    // it either subdivides further or stops subdividing.
    void subdivideSpaceAux(NodePtr ptr, const KdTree3& tree){

        real diagonal = halfPower[ptr->getDepth()];

        // We can assume that if we are at depth > 1,
        // that the parent of this node contained at least N_MIN points,
        // otherwise we would not have recursed.
        // So it is sufficient to search the paren't nodes for points 
        // contained inside this node.
        std::vector<Point> points;
        tree.rangeQuery(ptr->getCentroid(), ptr->getRadius(), points);

        // The empty ball case, no further subdivision
        if(points.size() == 0){
            
            
            return;
        }
        
        // Esnures there are enough points in the sphere
        while(points.size() < N_MIN){
            points.clear();
            ptr->incrementRadius();
            tree.rangeQuery(ptr->getCentroid(), ptr->getRadius(), points);
        }
        
        
        
    }


public:

    Octree() {
        // Root
        root = new Node(0, Vector3::ORIGIN);
    }



};

#endif