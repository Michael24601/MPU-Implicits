
#ifndef MATRIX_H
#define MATRIX_H

#include "precision.hpp"
#include <vector>


// A very simple matrix class, with limited operations.
// It is intended only for use in a Conjugate Gradient solver,
// for which it is well suited.
class Matrix{

private:

    std::vector<std::vector<real>> data;

public:

    Matrix(const std::vector<std::vector<real>>& data){
        assert((data.size() != 0 && data.size() == data[0].size())
            && "Matrix is not square\n");
        this->data = data;
    }


    Matrix(int n, real value = 0.0){
        data = std::vector<std::vector<real>>(n, std::vector<real>(n, value));
    }


    // Creates a matrix that is the outer product of the given vector
    Matrix(const std::vector<real>& v){

        assert(v.size() != 0 &&"Input can't be an empty vector\n");
        
        this->data.resize(v.size()); 
        for(int i = 0; i < v.size(); i++){
            this->data[i].resize(v.size());
            for(int j = 0; j < v.size(); j++){
                this->data[i][j] = v[i] * v[j];
            }
        }
    }


    int rows() const{
        return data.size();
    }


    int columns() const{
        return data[0].size();
    }


    // Scalar multiplication
    Matrix operator*(real s) const{
        Matrix result(*this);
        for(int i = 0; i < data.size(); i++){
            for(int j = 0; j < data[i].size(); j++){
                result.data[i][j] *= s;
            }
        }
        return result;
    }


    // Element-wise addition
    Matrix operator+(const Matrix& m) const{

        assert((m.data.size() != 0 && m.data.size() == data.size() 
            && m.data[0].size() == data[0].size()) 
            && "Matrix dimensions must match\n");  

        Matrix result(*this);
        for(int i = 0; i < data.size(); i++){
            for(int j = 0; j < data[i].size(); j++){
                result.data[i][j] += m.data[i][j];
            }
        }
        return result;
    }


    // Element-wise subtraction
    Matrix operator-(const Matrix& m) const{
        
        assert((m.data.size() != 0 && m.data.size() == data.size() 
            && m.data[0].size() == data[0].size())
            &&"Matrix dimensions must match\n");

        Matrix result(*this);
        for(int i = 0; i < data.size(); i++){
            for(int j = 0; j < data[i].size(); j++){
                result.data[i][j] -= m.data[i][j];
            }
        }
        return result;
    }


    // Vector multiplication
    std::vector<real> operator*(const std::vector<real>& v) const{
        
        assert((v.size() != 0 && v.size() == data[0].size())
            &&"Vector dimension must match\n");

        std::vector<real> result;
        result.resize(data.size());

        for(int i = 0; i < data.size(); i++){
            result[i] = 0;
            for(int j = 0; j < data[i].size(); j++){
                result[i] += data[i][j] * v[j];
            }
        }
        return result;
    }

};


#endif