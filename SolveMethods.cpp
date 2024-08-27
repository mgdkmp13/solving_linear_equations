#include "SolveMethods.h"
#include <ctime>
const int MAX_ITER = 1000;
const double BASE_RES_NORM = 2000000;

SolveMethods::SolveMethods() {
	this->A = Matrix();
	this->b = Matrix();
	this->resNorm = BASE_RES_NORM;
	this->time = 0;
	this->resNormVector = vector<double>();
	this->iterations = 0;
}

SolveMethods::SolveMethods(Matrix A, Matrix b) {
	this->A = A;
	this->b = b;
	this->resNorm = BASE_RES_NORM;
	this->time = 0;
	this->resNormVector = vector<double>();
	this->iterations = 0;
}

void SolveMethods::JacobiMethod() {
	this->resNorm = BASE_RES_NORM;
	this->time = 0;
	this->iterations = 0;
	this->resNormVector = vector<double>();
	clock_t start = clock();
	Matrix x = Matrix(this->b.rows, 1, 1);
	Matrix xNew = Matrix(this->b.rows, 1, 1);
	while (this->iterations < MAX_ITER && this->resNorm > pow(10, -9)) {
		for (int i = 0; i < this->A.rows; i++) {
			double sumHelp1 = 0.0;
			double sumHelp2 = 0.0;
			for (int j = 0; j < i; j++) {
				sumHelp1 += this->A.matrix[i][j] * x.matrix[j][0];
			}
			for (int j = i+1; j < this->A.rows; j++) {
				sumHelp2 += this->A.matrix[i][j] * x.matrix[j][0];
			}
			xNew.matrix[i][0] = (this->b.matrix[i][0] - sumHelp1 - sumHelp2) / this->A.matrix[i][i];
 		}
		this->calculateResNorm(xNew);
		x = xNew;
		this->iterations++;
	}
	clock_t end = clock();
	this->time = double(end - start) / CLOCKS_PER_SEC;
}

void SolveMethods::GaussSeidelMethod() {
	this->resNorm = BASE_RES_NORM;
	this->time = 0;
	this->iterations = 0;
	this->resNormVector = vector<double>();
	clock_t start = clock();
	Matrix x = Matrix(this->b.rows, 1, 1);
	Matrix xNew = Matrix(this->b.rows, 1, 1);
	while (this->iterations < MAX_ITER && this->resNorm > pow(10, -9)) {
		for (int i = 0; i < this->A.rows; i++) {
			double sumHelp1 = 0.0;
			double sumHelp2 = 0.0;
			for (int j = 0; j < i; j++) {
				sumHelp1 += this->A.matrix[i][j] * xNew.matrix[j][0];
			}
			for (int j = i + 1; j < this->A.rows; j++) {
				sumHelp2 += this->A.matrix[i][j] * x.matrix[j][0];
			}
			xNew.matrix[i][0] = (this->b.matrix[i][0] - sumHelp1 - sumHelp2) / this->A.matrix[i][i];
		}
		this->calculateResNorm(xNew);
		x = xNew;
		this->iterations++;
	}
	clock_t end = clock();
	this->time = double(end - start) / CLOCKS_PER_SEC;
}

void SolveMethods::LUMethod() {
	this->resNorm = BASE_RES_NORM;
	this->time = 0;
	this->resNormVector = vector<double>();
	clock_t start = clock();
	Matrix U = Matrix(this->A);
	Matrix L = Matrix(this->A.rows, this->A.cols);
	L.initializeEye();
	for (int i = 1; i < this->A.rows; i++) {
		for (int j = 0; j < i; j++) {
			L.matrix[i][j] = U.matrix[i][j] / U.matrix[j][j];
			for (int k = 0; k < this->A.rows; k++) {
				U.matrix[i][k] -= L.matrix[i][j] * U.matrix[j][k];
			}
		}
	}
	Matrix y = this->forwardSubstitution(L, this->b);
	Matrix x = this->backwardSubstitution(U, y);
	this->calculateResNorm(x);
	clock_t end = clock();
	this->time = double(end - start) / CLOCKS_PER_SEC;
}

Matrix SolveMethods::forwardSubstitution(Matrix mat, Matrix vec) {
	Matrix result = Matrix(mat.rows, vec.cols);
	for (int i = 0; i < mat.rows; i++) {
		double helpSum = 0.0;
		for (int j = 0; j < i; j++) {
			helpSum += mat.matrix[i][j] * result.matrix[j][0];
		}
		result.matrix[i][0] = (vec.matrix[i][0] - helpSum) / mat.matrix[i][i];
	}
	return result;
}

Matrix SolveMethods::backwardSubstitution(Matrix mat, Matrix vec) {
	Matrix result = Matrix(mat.rows, vec.cols);
	for (int i = mat.rows - 1; i >= 0 ; i--) {
		double helpSum = 0.0;
		for (int j = i + 1; j < mat.rows ; j++) {
			helpSum += mat.matrix[i][j] * result.matrix[j][0];
		}
		result.matrix[i][0] = (vec.matrix[i][0] - helpSum) / mat.matrix[i][i];
	}
	return result;
}

void SolveMethods::calculateResNorm(Matrix x) {
	Matrix helpMat = this->A.mulMat(x);
	helpMat = helpMat.subMat(this->b);
	this->resNorm = this->calculateNorm(helpMat);
	if (this->resNorm != INFINITY) {
		this->resNormVector.push_back(this->resNorm);
	}
	else {
		//this->resNorm = this->resNormVector.back();
		this->resNormVector.push_back(this->resNorm);
	}
}

double SolveMethods::calculateNorm(Matrix x) {
	double result = 0.0;
	for (int i = 0; i < x.rows; i++) {
		result += pow(x.matrix[i][0], 2);
	}
	return sqrt(result);
}

SolveMethods::~SolveMethods() {
}