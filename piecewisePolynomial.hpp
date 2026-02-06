
#ifndef PIECEWISE_POLYNOMIAL_H
#define PIECEWISE_POLYNOMIAL_H

#include "bivariateQuadratic.hpp"
#include "generalQuadric.hpp"
#include <optional>

class PiecewisePolynomial: public LocalFitFunction{

private:

    // One function for each cluster
    std::vector<BivariateQuadratic> functions;
    std::optional<GeneralQuadric> quadric;

    
    static void minimumDotProduct(const std::vector<Point>& points,
        int& n1, int& n2, real& min) {

        assert(points.size() >= 2 && "Points array is too small\n");
            
        n1 = 0;
        n2 = 1;
        min = points[n1].getNormal() * points[n2].getNormal();
        for(int i = 0; i < points.size(); i++){
            for(int j = i+1; j < points.size(); j++){
                real dot = points[i].getNormal() * points[j].getNormal();
                if(dot < min){
                    min = dot;
                    n1 = i;
                    n2 = j;
                }
            }
        }
    }


    // Given a cluster of points, divides them into two clusters
    // using the dot product test with two normals n1 and n2.
    static void clusterPoints(const std::vector<Point>& points,
        std::vector<Point>& p1, std::vector<Point>& p2, int n1, int n2){

        assert(points.size() >= 2 && "Points array is too small\n");
        
        p1.clear();
        p2.clear();

        for(int i = 0; i < points.size(); i++){
            real dot1 = points[i].getNormal() * points[n1].getNormal();
            real dot2 = points[i].getNormal() * points[n2].getNormal();
            if(dot1 >= dot2){
                p1.push_back(points[i]);
            }
            else{
                p2.push_back(points[i]);
            }
        }
    }


    // Calculates mean normal of cluster
    static Vector3 clusterMeanNormal(const std::vector<Point>& cluster){
        Vector3 mean(0, 0, 0);
        for(int i = 0; i < cluster.size(); i++){
            mean = mean + cluster[i].getNormal();
        }
        mean.normalize();
        return mean;
    }


    static Vector3 clusterMeanCentroid(const std::vector<Point>& cluster){
        Vector3 mean(0, 0, 0);
        for(int i = 0; i < cluster.size(); i++){
            mean = mean + cluster[i].getPoint();
        }
        mean = mean * (1.0 / cluster.size());
        return mean;
    }


    // Fits a bivariate quadratic on a cluster and pushes it into
    // the functions array.
    void addQuadratic(const std::vector<Point>& cluster, 
        const Vector3& centroid, const Vector3& normal, real radius){

        // We need at least 6 points to be able to fit the cluster
        BivariateQuadratic q(centroid, normal);
        if(cluster.size() >= 6){
            q.fit(cluster, centroid, radius);
        }
        functions.push_back(q);
    }


public:

    PiecewisePolynomial(){}

    bool fit(const std::vector<Point>& points, 
        const Vector3& centroid, real radius, 
        const Vector3& meanNormal, const std::vector<Vector3>& aux,
        const KdTree3& tree) {

        // We need at least 2 points
        assert(points.size() >= 2 && "Points array is too small\n");

        // Clears both
        functions.clear();
        quadric.reset();

        // First we must calculate the minimum dot product
        // between the normals of any two points.
        // In practice we only have 15-30 points, so I am sticking with 
        // the naive approach.
        int n1, n2;
        real min;
        minimumDotProduct(points, n1, n2, min);

        // If this min is >= 0.9, then we have one cluster
        if(min >= 0.9){
            addQuadratic(points, centroid, meanNormal, radius);
            return true;
        }

        // We can then compute the max absolute dot product with
        // the cross product of n1 and n2;
        Vector3 cross = (points[n1].getNormal() 
            % points[n2].getNormal()).normalized();

        real max = std::abs(points[0].getNormal() * cross);
        for(int i = 1; i < points.size(); i++){
            real dot = std::abs(points[i].getNormal() * cross);
            if(dot > max){
                max = dot;
            }
        }

        // If this max <= 0.7, then we have an edge (2 clusters)
        // and we use the normals from before n1 and n2 to build
        // the local space.
        if(max <= 0.7){
            // We then subdivide the points into two clusters,
            // each closer to n1 or n2.
            std::vector<Point> p1, p2;
            clusterPoints(points, p1, p2, n1, n2);

            // Each cluster must have at least 6 points, or we
            // keep the coefficients 0 in the function.
            addQuadratic(p1, clusterMeanCentroid(p1), clusterMeanNormal(p1), radius);
            addQuadratic(p2, clusterMeanCentroid(p2), clusterMeanNormal(p2), radius);

            return true;
        }

        // Otherwise, we have three clusters, and need to to cluster
        // them into three sets (see report for more details).
        std::vector<Point> p1, p2, p3;
        for(int i = 0; i < points.size(); i++){
            real dot1 = points[i].getNormal() * points[n1].getNormal();
            real dot2 = points[i].getNormal() * points[n2].getNormal();
            real dot3 = points[i].getNormal() * cross;

            if(dot1 < dot3 && dot2 < dot3){
                p3.push_back(points[i]);
            }
            else{
                if(dot1 >= dot2){
                    p1.push_back(points[i]);
                }
                else{
                    p2.push_back(points[i]);
                }
            }
        }

        // Now we repeat the minimum dot product test on the third cluster.
        // Unless of course p3 has less than two elements, meaning
        // we can't do that, and instead just go directly to the corner
        // case.
        int n31, n32;
        real min3;
        if(p3.size() >= 2){
            minimumDotProduct(p3, n31, n32, min3);
        }

        // If the test doesn't pass, we just have 3 clusters
        if(p3.size() < 2 || min3 >= 0.9){
            addQuadratic(p1, clusterMeanCentroid(p1), clusterMeanNormal(p1), radius);
            addQuadratic(p2, clusterMeanCentroid(p2), clusterMeanNormal(p2), radius);
            addQuadratic(p3, clusterMeanCentroid(p3), clusterMeanNormal(p3), radius);

            return true;
        }

        // Otherwise we subdivide the third cluster to get a degree 4
        // corner. First we need to check if the current surface
        // is even convex or concave, using the mean normal.
        bool isConvexConcave = true;
        for(int i = 0; i < p3.size(); i++){
            // If one of the normals is on the other side of the mean,
            // that means it is not convex or concave.
            if(meanNormal * p3[i].getNormal() < 0){
                isConvexConcave = false;
                break;
            }
        }

        // If yes, then we can fit 4 quadratics
        if(isConvexConcave){
            // We subdivide the third cluster
            std::vector<Point> p31, p32;
            clusterPoints(p3, p31, p32, n31, n32);

            addQuadratic(p1, clusterMeanCentroid(p1), clusterMeanNormal(p1), radius);
            addQuadratic(p2, clusterMeanCentroid(p2), clusterMeanNormal(p2), radius);
            addQuadratic(p31, clusterMeanCentroid(p31), clusterMeanNormal(p31), radius);
            addQuadratic(p32, clusterMeanCentroid(p32), clusterMeanNormal(p32), radius);

            return true;
        }

        // Otherwise, we fit a quadric
        quadric = GeneralQuadric();
        // Quadric fitting can always fail
        if(quadric->fit(points, aux, centroid, radius, tree)){
            return true;
        }
        else{
            return false;
        }

    }


    real evaluate(const Vector3& input) const override{

        assert(((functions.size() > 0 && functions.size() <= 4) 
            || quadric) && "The number of functions is not correct\n");

        // We just need to use max
        if(functions.size() == 1){
            return functions[0].evaluate(input);
        }
        else if (functions.size() == 2){
            return std::max(functions[0].evaluate(input),
                functions[1].evaluate(input));
        }
        else if(functions.size() == 3) {
            return std::max(functions[0].evaluate(input),
                std::max(functions[1].evaluate(input), 
                    functions[2].evaluate(input)));
        }
        else if(functions.size() == 4) {
            return std::max(functions[0].evaluate(input),
                std::max(functions[1].evaluate(input), 
                    std::max(functions[2].evaluate(input), 
                        functions[3].evaluate(input))));
        }

        return quadric->evaluate(input);

    }

    
    Vector3 evaluateGradient(const Vector3& input) const override{

        assert(((functions.size() > 0 && functions.size() <= 4) 
            || quadric) && "The number of functions is not correct\n");

        // We just need to use min() (since we have 
        // a function that is negative on the inside of the surface
        // and poisitive outside).
        if(functions.size() == 1){
            return functions[0].evaluateGradient(input);
        }
        else if (functions.size() == 2){
            real q1 = functions[0].evaluate(input);
            real q2 = functions[1].evaluate(input);

            if(q1 >= q2){
                return functions[0].evaluateGradient(input);
            }
            else{
                return functions[1].evaluateGradient(input);
            }
        }
        else if(functions.size() == 3) {
            real q1 = functions[0].evaluate(input);
            real q2 = functions[1].evaluate(input);
            real q3 = functions[2].evaluate(input);

            if(q1 >= q2 && q1 >= q3){
                return functions[0].evaluateGradient(input);
            }
            else if (q2 >= q3 && q2 > q1){
                return functions[1].evaluateGradient(input);
            }
            else{
                return functions[2].evaluateGradient(input);
            }
        }
        else if(functions.size() == 4) {
            real q1 = functions[0].evaluate(input);
            real q2 = functions[1].evaluate(input);
            real q3 = functions[2].evaluate(input);
            real q4 = functions[3].evaluate(input);

            if(q1 >= q2 && q1 >= q3 && q1 >= q4){
                return functions[0].evaluateGradient(input);
            }
            else if (q2 >= q3 && q2 >= q4 && q2 > q1){
                return functions[1].evaluateGradient(input);
            }
            else if (q3 >= q4 && q3 > q2 && q3 > q1){
                return functions[2].evaluateGradient(input);
            }
            else{
                return functions[3].evaluateGradient(input);
            }
        }

        return quadric->evaluateGradient(input);
    }

};

#endif