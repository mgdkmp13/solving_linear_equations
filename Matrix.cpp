#include "Matrix.h"
const int BASE_SIZE = 2;

Matrix::Matrix() {
	for (int i = 0; i < BASE_SIZE; i++) {
		vector<double> newV(BASE_SIZE, 0);
		this->matrix.push_back(newV);
	}
	this->rows = BASE_SIZE;
	this->cols = BASE_SIZE;
}

Matrix::Matrix(int m, int n) {
	for (int i = 0; i < m; i++) {
		vector<double> newV(n, 0);
		this->matrix.push_back(newV);
	}
	this->rows = m;
	this->cols = n;
}

Matrix::Matrix(int m, int n, double base) {
	for (int i = 0; i < m; i++) {
		vector<double> newV(n, base);
		this->matrix.push_back(newV);
	}
	this->rows = m;
	this->cols = n;
}

Matrix::Matrix(const Matrix& other) : matrix(other.matrix), rows(other.rows), cols(other.cols) {}

void Matrix::initialize(double a1, double a2, double a3) {
	for (int y = 0; y < this->rows; y++) {
		for (int x = 0; x < this->cols; x++) {
			if (x == y) {
				this->matrix[y][x] = a1;
			}
			else if (abs(y-x)==1) {
				this->matrix[y][x] = a2;
			}
			else if (abs(y-x)==2) {
				this->matrix[y][x] = a3;
			}
		}
	}
}

void Matrix::initializeBVector(double F) {
	for (int i = 0; i < this->rows; i++) {
		this->matrix[i][0] = sin(i * (F + 1.0));
	}
}

void Matrix::initializeEye() {
	for (int i = 0; i < this->rows; i++) {
		this->matrix[i][i] = 1.0;
	}
}

Matrix Matrix::mulMat(Matrix mat) {
	int help = 0;
	Matrix resultMat(this->rows, mat.cols);
	for (int y = 0; y < this->rows; y++) {
		for (int x = 0; x < mat.cols; x++) {
			for (int i = 0; i < this->cols; i++) {
				resultMat.matrix[y][x] += this->matrix[y][i] * mat.matrix[i][x];
			}
		}
	}
	return resultMat;
}

Matrix Matrix::subMat(Matrix mat) {
	Matrix resultMat(this->rows, this->cols, 0);
	for (int y = 0; y < this->rows; y++) {
		for (int x = 0; x < this->cols; x++) {
			resultMat.matrix[y][x] = this->matrix[y][x] - mat.matrix[y][x];
		}
	}
	return resultMat;
}


void Matrix::print() {
	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			if (this->matrix[i][j] < 0) {
				cout << this->matrix[i][j] << " ";
			}
			else {
				cout << " " << this->matrix[i][j] << " ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

Matrix::~Matrix() {
}