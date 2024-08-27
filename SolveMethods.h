#pragma once
#include <iostream>
#include <vector>
#include "Matrix.h"
using namespace std;


class SolveMethods {
public:
    Matrix A;
    Matrix b;
    double resNorm;
    double time;
    short iterations;
    vector<double> resNormVector;
    SolveMethods();
    SolveMethods(Matrix A, Matrix b);
    double calculateNorm(Matrix x);
    void calculateResNorm(Matrix x);
    void JacobiMethod();
    void GaussSeidelMethod();
    void LUMethod();
    Matrix forwardSubstitution(Matrix mat, Matrix vec);
    Matrix backwardSubstitution(Matrix mat, Matrix vec);
    ~SolveMethods();
};
