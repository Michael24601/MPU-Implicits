
#ifndef SOLVER_TEST_H
#define SOLVER_TEST_H

#include "../../conjugateGradient.hpp"
#include <random>
#include <iostream>

// This function generates a random SPD matrix with dimension n
Matrix generateSPDMatrix(int n, int seed){

    // Initially filled with 0s
    std::vector<std::vector<real>> data(n, std::vector<real>(n, 0.0));

    // Random value in between 0.5, and -0.5.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> distribution(-0.5, 0.5);

    // We then just fill the upper triangle and then mirror the
    // result on the lower triangle.
    for (int i = 0; i < n; i++){
        real rowSum = 0.0;

        for (int j = i+1; j < n; j++){
            data[i][j] = data[j][i] = distribution(rng);

            rowSum += std::abs(data[i][j]);
        }

        // Then we just ensure that the diagonals are dominant
        data[i][i] = rowSum + 0.5 + std::abs(distribution(rng));
    }

    Matrix result(data);
    return result;
}


// Generates the true vector of unknowns
std::vector<real> generateVector(int n, real min, real max, int seed) {
    std::vector<real> result(n);

    // Random value
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> distribution(min, max);

    for (int i = 0; i < n; i++){
        result[i] = distribution(rng);
    }
    return result;
}


void solverTest(){
    // We generate A and c
    Matrix A = generateSPDMatrix(20, 123456);
    std::vector<real> c = generateVector(20, -1.0, 1.0, 456789);

    // We then multiply them together to get b
    std::vector<real> b = A * c;

    // We then try to get c back using CG
    std::vector<real> cgOutput = ConjugateGradient::solve(A, b);
    
    // We can print Ax - b (normalized by b), which should be small
    real norm = sqrtReal(Matrix::dot(b, b));
    std::vector<real> error = Matrix::subtract(A*cgOutput, b);
    real residualNorm = sqrtReal(Matrix::dot(error, error));
    residualNorm /= norm;
    std::cout << residualNorm << "\n";
}

#endif