
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
#include "bivariateQuadratic.hpp"
#include "kdTree3.hpp"
#include "cubeOffset.hpp"
#include "weightFunction.hpp"

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

        
        // Calculates the maximum angle theta between
        // the normals of the points and the normalized arithmetic
        // mean of the points.
        void computeMaxTheta(const std::vector<Point>& points, 
            real& cosTheta, Vector3& meanNormal) const {

            if(points.empty()){
                cosTheta = 0.0;
                meanNormal = Vector3::ORIGIN;
                return;
            }

            // We need to calculate the weighted average,
            // using the B-spline weights.
            Vector3 mean = Vector3::ORIGIN;
            for(int i = 0 ; i < points.size(); i++){
                mean = mean + points[i].getNormal() 
                    * evaluateWeightFunction(points[i].getPoint());
            }
            // It's enough to normalized to average out (since we
            // need to normalize anyway)
            mean.normalize();

            // Max angle is the minimum dot product (no need
            // to find the angle, using cos(theta) gives us as much
            // information since the angle is restricted to [0, PI]).

            real minDotProduct = 9999.0;
            for(int i = 0; i < points.size(); i++){
                minDotProduct = std::min(minDotProduct, 
                    mean * points[i].getNormal());
            }
            cosTheta = minDotProduct;
            meanNormal = mean;
        } 


    public:

        Node(int depth, const Vector3& centroid) : 
            depth{depth}, centroid{centroid}, leaf{true},
            child{nullptr, nullptr, nullptr, nullptr, nullptr, 
                nullptr, nullptr, nullptr}, localFitFunction{nullptr}{

            // By default the radius is ALPHA * diagonal()    
            radius = ALPHA * getDiagonal();
         }

        
        ~Node(){
            deallocateLocalFunction();
        }


        void deallocateLocalFunction(){
            if(localFitFunction){
                delete localFitFunction;
                localFitFunction = nullptr;
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
            radius = radius + LAMBDA * (ALPHA * getDiagonal());
        }


        bool isLeaf() const{
            return leaf;
        }


        void setRadius(real radius){
            this->radius = radius;
        }


        // The diagonal of each subdivision is some power (1/2)^n.
        real getDiagonal() const{
            return halfPower[depth];
        }


        // Subdivides at this node by initializing 8 children
        bool subdivide(){
            // If the current node can't be subdivided for any reason,
            // we return false
            if(depth == MAX_DEPTH || !leaf){
                return false;
            }

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
            
            return true;
        }


        // Generates the auxiliary points of a node (corners and centroid).
        std::vector<Vector3> generateAuxiliaryPoints() const{
            // The auxiliary points are the corners and the centroid
            std::vector<Vector3> auxiliary(9);
            auxiliary[0] = centroid;
            // We need half of this cube's diagonal
            real offset = halfPower[depth+1] * INV_SQRT3;
            
            for(int i = 0; i < 8; i++){
                auxiliary[i+1] = centroid + (cubeCenterOffset[i] * offset);
            }
            return auxiliary;
        }

    
        // Assuming that this node contains the given points,
        // this function fits a local function to the points.
        // Note that this function assumes that points is not empty.
        // This function can return false if no fitting was able to be
        // done.
        bool fitLocalFunction(const std::vector<Point>& points, 
            const KdTree3& tree){

            // Only fitted on leafs
            if(!leaf){
                return false;
            }

            // This just computes the maximum angle theta
            // between the normals and the mean normal,
            // as well as the mean normal;
            real cosTheta;
            Vector3 meanNormal;
            computeMaxTheta(points, cosTheta, meanNormal);

            // If cosTheta < 0, then theta > PI/2 (at least if we 
            // restrict theta to [0, PI])
            if(cosTheta <= 0){

                // General quadric case
                GeneralQuadric* quadric = new GeneralQuadric();
                
                // The auxiliary points are the corners and the centroid
                std::vector<Vector3> auxiliary = generateAuxiliaryPoints();

                if(quadric->fit(points, auxiliary, centroid, radius, tree)){
                    this->localFitFunction = quadric;
                    return true;
                }
                else{
                    delete quadric;
                    return false;
                }
            }
            else if(true){
                BivariateQuadratic* quadratic = 
                    new BivariateQuadratic(centroid, meanNormal);
                quadratic->fit(points, centroid, radius);
                this->localFitFunction = quadratic;
                return true;
            }
            else{
                // Edges and corners
                return false;
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
                std::cerr << "Local function undefined on node\n";
                return 0.0;
            }
        }


        real evaluateErrorApproximation(
            const std::vector<Point>& points) const{

            // In general only leafs should have local functions fitted,
            // so we need to ensure we only call this on leafs.
            if(localFitFunction){
                // We scale the error by teh diagonal
                return localFitFunction->approximationError(points);
            }
            else{
                throw std::runtime_error("Local function undefined on node\n");
            }
        }


        // Evaluates the b-spline of the local function at the
        // given point.
        real evaluateWeightFunction(const Vector3& p) const {
            // Just calls the B-spline evaluate function
            return WeightFunction::evaluate(centroid, radius, p);
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

        // In the empty ball case, no further subdivision is done.
        // We increase the radius as usual, but we don't subdivide it
        // later.
        bool emptyFlag = false;
        if(points.size() == 0){
            emptyFlag = true;
        }

        // If we have some points in the sphere, we ensure we have
        // enough points, by incrementing the radius.
        while(points.size() < N_MIN){
            points.clear();
            ptr->incrementRadius();
            tree.rangeQuery(ptr->getCentroid(), ptr->getRadius(), points);
        }

        // Once we have enough points, we fit a local function.
        // If the function returns false, it means it tried to fit
        // a local function but for some reason the result is not
        // perfectly usable.
        // We thus prioritize subdividing (regardless of error
        // approximation), unless the MAX_DEPTH has been reached
        // or the sphere was originally empty, in which case we
        // don't subdivide.
        if(!ptr->fitLocalFunction(points, tree)){

            // We immediatly subdivide if we can subdivide
            if(!emptyFlag && ptr->subdivide()){
                points.clear();
                for(int i = 0; i < 8; i++){
                    // Then repeat on the children
                    subdivideSpaceAux(ptr->getChild(i), tree);
                }
                return;
            }
 
            // Otherwise, if we can't subdivide because we reached
            // the max depth for example, the function remains empty.
            // This will cause a runtime error if we try to evaluate the
            // leaf node.

            return;
        }

        // If this flag is true, no subdivision is needed, so no need
        // for an error approximation.
        if(emptyFlag){
            return;
        }

        // If this error metric is too large, we have to disregard
        // the local fit function and subdivide.
        real error = ptr->evaluateErrorApproximation(points);
        if(error > EPSILON_ZERO){

            // If the error is too large, we delete the function and
            // subdivide, unless we can't (if MAX_DEPTH is reached)
            // then we just keep the function as is.
            if(ptr->subdivide()){
                points.clear();
                ptr->deallocateLocalFunction();
                // Then repeat on the children
                for(int i = 0; i < 8; i++){
                    subdivideSpaceAux(ptr->getChild(i), tree);
                }
            }
            return;
        }
        
        // Otherwise we keep the function and the node remains a leaf
    }


    void evaluateAux(const Vector3& p, NodePtr ptr, real& sum, 
        real& normalizationFactor){
        
        if(!ptr){
            return;
        }
        
        if((p - ptr->getCentroid()).lengthSquared() <= 
            ptr->getRadius() * ptr->getRadius()){

            // If it is a leaf, we evaluate it
            if(ptr->isLeaf()){

                real spline = ptr->evaluateWeightFunction(p);
                sum += ptr->evaluateLocalFitFunction(p) * spline;
                normalizationFactor += spline;
            }
            // Otherwise we check its children
            else{
                // Then we visit the children
                for(int i = 0; i < 8; i++){
                    evaluateAux(p, ptr->getChild(i), sum, normalizationFactor);
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
        // In case the point is outside the root cube
        if(std::abs(p.x()) > INV_SQRT3 * 0.5 
            || std::abs(p.y()) > INV_SQRT3 * 0.5
            || std::abs(p.z()) > INV_SQRT3 * 0.5){
            return 0.0;
        }
        
        // Sum and normalization factor
        real sum = 0.0;
        real normalizationFactor = 0.0;
        // We will go through the tree, exploring branches only
        // if the sphere contains the given point.
        // When we get to a leaf, and it contains the point,
        // we evaluate the weight and the local function associated
        // with it, and add it to the sum.
        // We also keep track of the sum of weight functions (those
        // that are non-zero) so we can normalize the result.
        evaluateAux(p, root, sum, normalizationFactor);
        if(normalizationFactor > 0.0){
            sum /= normalizationFactor;
        }
        else{
            sum = 0.0;
        }
        return sum;
    }

};

#endif