
#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include "precision.hpp"
#include "vector3.hpp"
#include "localFitFunction.hpp"
#include "generalQuadric.hpp"
#include "kdTree3.hpp"
#include "cubeOffset.hpp"


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
        // The local fit function, should only be fitted
        // on a leaf node. If a node is subdivided we should delete
        // the local fit function.
        LocalFitFunction* localFitFunction;
        
        // The 8 children
        Node* child[8] = {nullptr};
        bool leaf;

    public:

        Node(int depth, const Vector3& centroid) : 
            depth{depth}, centroid{centroid}, leaf{true},
            child{nullptr, nullptr, nullptr, nullptr, nullptr, 
                nullptr, nullptr, nullptr}{

            // By default the radius is ALPHA * diagonal()    
            radius = ALPHA * getDiagonal();
         }

        
        ~Node(){
            deallocateLocalFunction();
        }


        void deallocateLocalFunction(){
            if(localFitFunction){
                delete localFitFunction;
            }
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


        Node* getChild(int index){
            if(index < 0 || index >= 8){
                throw std::invalid_argument("Child index does not exist\n");
            }
            return child[index];
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
                // Ceases to be a leaf
                leaf = false;

                // The diagonal for this cube is halfPower[depth],
                // so the offset diagonal is at the quarter, that is is,
                // halfPower[depth+2], which we know exists since we
                // haven't reached the max depth.
                // We multiply it by 1/sqrt3 to get the offset in x,
                // y and z.
                real offset = halfPower[depth+2] * INV_SQRT3;

                for(int i = 0; i < 8; i++){
                    child[i] = new Node(depth+1, 
                        centroid + (cubeCenterOffset[i] * offset));
                }
            }
        }

    
        // Assuming that this node contains the given points,
        // this function fits a local function to the points.
        // Note that this function assumes that points is not empty.
        void fitLocalFunction(const std::vector<Point>& points){

            // Only fitted on leafs
            if(!leaf){
                return;
            }

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
                    corners[i] = centroid + (cubeCenterOffset[i] * offset);
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

        // Evaluates the local function at the given point
        real evaluateLocalFitFunction(const Vector3& p) const{
            // In general only leafs should have local functions fitted,
            // so we need to ensure we only call this on leafs.
            if(localFitFunction){
                return localFitFunction->evaluate(p);
            }
            else{
                throw std::exception("Local function undefined on node\n");
            }
        }


        real evaluateErrorApproximation(
            const std::vector<Point>& points) const{

            // In general only leafs should have local functions fitted,
            // so we need to ensure we only call this on leafs.
            if(localFitFunction){
                return localFitFunction->approximationError(points);
            }
            else{
                throw std::exception("Local function undefined on node\n");
            }
        }


        // Evaluates the b-spline of the local function at the
        // given point
        real evaluateQuadraticBSpline(const Vector3& p) const {
            // TODO: finish this.
            return 1.0;
        }

    };


    typedef Node* NodePtr;
    NodePtr root;


    // Given a node representing a cube in the octree, this function
    // computes the local fit function, and if some criteria is met,
    // it either subdivides further or stops subdividing.
    void subdivideSpaceAux(NodePtr ptr, const KdTree3& tree){

        // As soon as a node is created, its radius is set to
        // ALPHA * diagonal (which is fixed given the depth), 
        // and the centroid is given in the constructor, and is
        // also fixed based on the paren't centroid and diagonal.
        // We can assume ptr has been initialized, if not,
        // we return.
        if(!ptr) return;

        // We assume a pointer that made it to this point is a leaf.
        // If not, then 

        // We check how many points are in the sphere by default
        std::vector<Point> points;
        tree.rangeQuery(ptr->getCentroid(), ptr->getRadius(), points);

        // In the empty ball case, no further subdivision i done,
        // we fit the local function and move on.
        if(points.size() == 0){
            ptr->fitLocalFunction(points);
            return;
        }
        
        // If we have some points in the sphere, we ensure we have
        // enough points, by incrementing the radius.
        while(points.size() < N_MIN){
            points.clear();
            ptr->incrementRadius();
            tree.rangeQuery(ptr->getCentroid(), ptr->getRadius(), points);
        }
        
        // Once we have enough points, we fit a local function.
        ptr->fitLocalFunction(points);
        
        // If this error metric is too large, we have to disregard
        // the local fit function and subdivide.
        if(ptr->evaluateErrorApproximation(points) > EPSILON_ZERO){
            ptr->deallocateLocalFunction();
            ptr->subdivide();
            for(int i = 0; i < 8; i++){
                subdivideSpaceAux(ptr->getChild(i), tree);
            }
        }        

        // Otherwise we keep the function and the node remains a leaf
    }


    void evaluateAux(const Vector3& p, NodePtr ptr, float& sum, 
        float& factor){
        
        if(!ptr){
            return;
        }
        
        if((p - ptr->getCentroid()).lengthSquared() <= 
            ptr->getRadius() * ptr->getRadius()){

            // If it is a leaf, we evaluate it
            if(ptr->isLeaf()){

                real spline = ptr->evaluateQuadraticBSpline(p);
                sum += ptr->evaluateLocalFitFunction(p) * spline;
                factor *= spline;
            }
            // Otherwise we check its children
            else{
                // Then we visit the children
                for(int i = 0; i < 8; i++){
                    evaluateAux(p, ptr->getChild(i), sum, factor);
                }
            }
        }
    }


public:

    Octree() {
        // Root
        root = new Node(0, Vector3::ORIGIN);
    }


    // Subdivides the space given points arranged in a kd-tree.
    // This assumes the points are in a cube centered at the origin
    // with a diagonal of size 1.
    void subdivideOctree(const KdTree3& tree){
        // We can't do this iteratively since we may need
        // to subdivide the nodes.
        subdivideSpaceAux(root, tree);
    }


    real evaluate(const Vector3& p){
        // Sum and normalization factor
        float sum = 0.0;
        float factor = 1.0;
        // We will go through the tree, exploring branches only
        // if the sphere contains the given point.
        // When we get to a leaf, and it contains the point,
        // we evaluate the weight and the local function associated
        // with it, and add it to the sum.
        // We also keep track of the sum of weight functions (those
        // that are non-zero) so we can normalize the result.
        evaluateAux(p, root, sum, factor);
        sum /= factor;
        return sum;
    }

};

#endif