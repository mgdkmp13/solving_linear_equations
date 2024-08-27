#pragma once
#include <iostream>
#include <vector>
using namespace std;


class Matrix {
public:
    vector<vector<double>> matrix;
    int rows;
    int cols;
    Matrix();
    Matrix(int m, int n);
    Matrix(int m, int n, double base);
    Matrix(const Matrix& other);
    void initialize(double a1, double a2, double a3);
    void initializeBVector(double F);
    void initializeEye();
    Matrix mulMat(Matrix mat);
    Matrix subMat(Matrix mat);
    void print();
    ~Matrix();
};