
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
#include "piecewisePolynomial.hpp"
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
        
        // This is the support of the entire subtree originating
        // at this node. It is a cuboid, so it will in practice
        // be larger than the actual support (union of spheres)
        std::pair<Vector3, Vector3> subtreeSupport;
        
        // The 8 children
        Node* child[8];
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

            assert((index >= 0 && index < 8)
                && "Child index does not exist\n");
                
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


        bool isFitted() const{
            return localFitFunction != nullptr;
        }


        // Returns a constant pointer (can't change what it
        // points to, can't change the object it points to)
        const LocalFitFunction* const getLocalFunction() const{
            return localFitFunction;
        }


        const std::pair<Vector3, Vector3>& getSubtreeSupport() const {
            return subtreeSupport;
        }


        // Sets the support that the entire subtree has originating at'
        // this current node. This support is a cuboid.
        // This is not the extact support, since the support will
        // be a union of spheres, but it cotnains the support and
        // prunes most points outside of it.
        // Note that it also sets the descendants' support as it works
        // from the bottom up.
        void computeSubtreeSupport() {

            std::pair<Vector3, Vector3> bounds(
                centroid - Vector3(radius),
                centroid + Vector3(radius)
            );
            
            if(!leaf){
                for(int i = 0; i < 8; i++){
                    child[i]->computeSubtreeSupport();
                    std::pair<Vector3, Vector3>& b = child[i]->subtreeSupport;
                    // We then possibly expand the bounds based on the
                    // children.
                    for(int c = 0; c < 3; c++){
                        bounds.first[c] = std::min(b.first[c], bounds.first[c]);
                        bounds.second[c] = std::max(b.second[c], bounds.second[c]);
                    }
                }
            }

            this->subtreeSupport = bounds;
        }


        // Subdivides at this node by initializing 8 children
        bool subdivide(){
            // If the current node can't be subdivided for any reason,
            // we return false
            if(depth == Octree::MAX_DEPTH || !leaf){
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
            // restrict theta to [0, PI]).
            if(cosTheta <= 0 && 
                (points.size() >= 2 * N_MIN || !USE_PIECEWISE_POLYNOMIALS)){

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
            else if(points.size() >= 2 * N_MIN || !USE_PIECEWISE_POLYNOMIALS){
                BivariateQuadratic* quadratic = 
                    new BivariateQuadratic(centroid, meanNormal);
                quadratic->fit(points, centroid, radius);
                this->localFitFunction = quadratic;
                return true;
            }
            else{
                // Edges and corners
                // The auxiliary points are the corners and the centroid
                std::vector<Vector3> auxiliary = generateAuxiliaryPoints();
                PiecewisePolynomial* p = new PiecewisePolynomial();
                if(p->fit(points, centroid, radius, meanNormal, 
                    auxiliary, tree)){
                    this->localFitFunction = p;
                    return true;
                }
                else{
                    delete p;
                    return false;
                }
            }
        }

        // Evaluates the local function at the given point
        real evaluateLocalFitFunction(const Vector3& p) const{
            // In general only leafs should have local functions fitted,
            // so we need to ensure we only call this on leafs.
            // However, it is possible that no function was able to be
            // fitted; this doesn't throw an error, but it does alert
            // the user, and returns 0.0.
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
            assert(localFitFunction != nullptr 
                && "Local function not defined");

            // We scale the error by teh diagonal
            return localFitFunction->approximationError(points);
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
            // The user will be alerted if an attempt is made to
            // evaluate it.

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

            // Note, we don't deallocate the parent function, since
            // we may need it as a fallback in case the children could not
            // be fitted.
            if(ptr->subdivide()){
                points.clear();
                // Then repeat on the children
                for(int i = 0; i < 8; i++){
                    subdivideSpaceAux(ptr->getChild(i), tree);
                }
            }
            return;
        }
        
        // Otherwise we keep the function and the node remains a leaf
    }


    void evaluateAux(const Vector3& p, NodePtr ptr, 
        const LocalFitFunction* const parentFunction, real& sum, 
        real& normalizationFactor){
        
        if(!ptr){
            return;
        }
        
        // If the node is a leaf and the node is in its sphere
        // support, then we can evalute it.
        if(ptr->isLeaf()){
            // If it is within the bounds, we evaluate it, otherwise
            // we return early.
            if((p - ptr->getCentroid()).lengthSquared() <= 
                ptr->getRadius() * ptr->getRadius()){

                real spline = ptr->evaluateWeightFunction(p);

                // If we can't evaluate the current local function,
                // because it was not fitted, we evaluate the parent,
                // or closest ancestor that exists.
                real eval;
                if(ptr->isFitted()){
                    eval = ptr->evaluateLocalFitFunction(p);
                }
                else{
                    if(parentFunction){
                        eval = parentFunction->evaluate(p);
                    }
                    else{
                        eval = 0.0;
                        std::cerr << "Could not evaluate the cell or "
                            << "any of its ancestors\n";
                    }
                }

                sum += eval * spline;
                normalizationFactor += spline;
            }
            return;
        }
        // If it is not a leaf, then we need to check its descendants,
        // but only if the descendant support contains the point.
        // Note that this function that checks for descendant support
        // overestimates the support, it may return true even though
        // the point is not in any of the descendants, since it is 
        // a cuboid, and the support is a union of spheres.
        // We can't just prune the subtree if the parent sphere doesn't
        // contain the point as the child spheres aren't necessarily contained
        // in the parent spheres.

        // We always keep track of the closest ancestor local function
        // (which may need to be evaluated in case the leaf doesn't
        // have a local function).
        else if(p.inBounds(ptr->getSubtreeSupport())){
            const LocalFitFunction* const f = (ptr->isFitted() ? 
                ptr->getLocalFunction(): parentFunction);
                
            for(int i = 0; i < 8; i++){
                evaluateAux(p, ptr->getChild(i), f, sum, normalizationFactor);
            }
        }
    }


    int getDepthAux(const Vector3& p, NodePtr ptr){
        if(!ptr || ptr->isLeaf()) {
            return 0;
        }

        // If the node is not a leaf, we need to know
        // which child p is in based on x, y, and z.
        int i = (p.x() > ptr->getCentroid().x() ? 1 : 0);
        int j = (p.y() > ptr->getCentroid().y() ? 1 : 0);
        int k = (p.z() > ptr->getCentroid().z() ? 1 : 0);
        return getDepthAux(p, ptr->getChild(inverseOffset[i][j][k])) + 1;
    }


public:

    // The minimum number of points in a leaf node
    static unsigned int N_MIN;

    // The factor used in the radius computation
    static real ALPHA;

    // The step size when incrementing the radius
    static real LAMBDA;

    // The maximum depth the octree can go
    static unsigned int MAX_DEPTH;

    // The maximum accepted approximation error
    // (Must be careful with this one, too small a value and the octree
    // never stops subsdividing and the max depth is always reached)
    static real EPSILON_ZERO;

    // Controls whether we use piecewise polynomials (sharp features)
    static bool USE_PIECEWISE_POLYNOMIALS;
    

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

        // Because the subtree supports never change, it is best
        // to precompute them here than each time we evaluate.
        root->computeSubtreeSupport();
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
        evaluateAux(p, root, nullptr, sum, normalizationFactor);
        if(normalizationFactor > 0.0){
            sum /= normalizationFactor;
        }
        else{
            sum = 0.0;
        }
        return sum;
    }


    // Returns the depth of subdivision of the given point.
    // This can be used to color the area, or to visualize the octree.
    int getDepth(const Vector3& p){
        return getDepthAux(p, root);
    }

};

#endif