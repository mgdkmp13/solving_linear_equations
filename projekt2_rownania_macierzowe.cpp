#include <iostream>
#include <vector>
#include <cmath>
#include <string>	
#include <fstream>
#include "Matrix.h"
#include "SolveMethods.h"
using namespace std;

const int C = 9;
const int D = 5;
const double E = 1.0;
const double F = 3.0;

void saveToCSV(short size, vector<double>& data, string& filename) {
	ofstream file(filename);

	if (!file.is_open()) {
		cout << "File opening error" << endl;
		return;
	}
	file << size << endl;
	file << data.size() << endl;
	for (size_t i = 0; i < data.size(); ++i) {
		file << data[i] << endl;
	}

	file.close();
}

void saveToCSV(double time, string& filename) {
	ofstream file(filename, ios::app);

	if (!file.is_open()) {
		cout << "File opening error" << endl;
		return;
	}
	file << time << endl;
	file.close();
}

void zadA_B(int N, Matrix b) {
	cout << endl << "%%%%%%%%%% ZAD B %%%%%%%%%%" << endl << endl;
	Matrix A(N, N);
	A.initialize(5.0 + E, -1.0, -1.0);
	SolveMethods solve(A, b);
	
	solve.JacobiMethod();
	cout << "~~ JACOBI ~~" << endl;
	cout << "Iterations: " << solve.iterations << endl;
	cout << "Residual norm: " << solve.resNorm << endl;
	cout << "Time: " << solve.time << "s" << endl << endl;
	string fileName = "jacobiA.txt";
	saveToCSV(solve.iterations, solve.resNormVector, fileName);
	
	solve.GaussSeidelMethod();
	cout << "~~ GAUSS-SEIDEL ~~" << endl;
	cout << "Iterations: " << solve.iterations << endl;
	cout << "Residual norm: " << solve.resNorm << endl;
	cout << "Time: " << solve.time << "s" << endl << endl;
	fileName = "gauss-seidelA.txt";
	saveToCSV(solve.iterations, solve.resNormVector, fileName);
}

void zadC_D(int N, Matrix b) {
	cout << endl << "%%%%%%%%%% ZAD C %%%%%%%%%%" << endl << endl;
	Matrix A(N, N);
	A.initialize(3.0, -1.0, -1.0);
	SolveMethods solve(A, b);
	
	solve.JacobiMethod();
	cout << "~~ JACOBI ~~" << endl;
	cout << "Iterations: " << solve.iterations << endl;
	cout << "Residual norm: " << solve.resNorm << endl;
	cout << "Time: " << solve.time << "s" << endl << endl;
	string fileName = "jacobiC.txt";
	saveToCSV(solve.iterations, solve.resNormVector, fileName);
	
	solve.GaussSeidelMethod();
	cout << "~~ GAUSS-SEIDEL ~~" << endl;
	cout << "Iterations: " << solve.iterations << endl;
	cout << "Residual norm: " << solve.resNorm << endl;
	cout << "Time: " << solve.time << "s" << endl<<endl;
	fileName = "gauss-seidelC.txt";
	saveToCSV(solve.iterations, solve.resNormVector, fileName);
	

	cout << endl << "%%%%%%%%%% ZAD D %%%%%%%%%%" << endl;
	solve.LUMethod();
	cout << "~~ LU ~~" << endl;
	cout << "Residual norm: " << solve.resNorm << endl;	
	cout << "Time: " << solve.time << "s" << endl;
	fileName = "LU_C.txt";
	saveToCSV(solve.iterations, solve.resNormVector, fileName);
}

void zadE() {
	cout << endl << "%%%%%%%%% ZAD E %%%%%%%%%" << endl;
	const int samples = 12;
	int sizes[samples] = { 100, 200, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };

	ofstream file("jacobiE.txt");
	file.close();
	file.open("gauss-seidelE.txt");
	file.close();
	file.open("LU_E.txt");
	file.close();

	for (int i = 0; i < samples; i++) {
		int N = sizes[i];
		cout << endl << "%~%~%~%~% N = " << N << " %~%~%~%~%" << endl << endl;
		Matrix A(N, N);
		A.initialize(5.0 + E, -1.0, -1.0);
		Matrix b(N, 1);
		b.initializeBVector(F);
		SolveMethods solve(A, b);
		
		solve.JacobiMethod();
		cout << "~~ JACOBI ~~";
		cout << endl << "Time: " << solve.time << "s" << endl;
		cout << "Residual norm:" << solve.resNorm << endl;
		string fileName = "jacobiE.txt";
		saveToCSV(solve.time, fileName);
		
		solve.GaussSeidelMethod();
		cout << endl << "~~ GAUSS-SEIDEL ~~";
		cout << endl << "Time: " << solve.time << "s" << endl;
		cout << "Residual norm:" << solve.resNorm << endl;
		fileName = "gauss-seidelE.txt";
		saveToCSV(solve.time, fileName);
		
		solve.LUMethod();
		cout << endl << "~~ LU ~~" << endl;
		cout << "Time: " << solve.time << "s" << endl;
		cout << "Residual norm:" << solve.resNorm << endl;
		fileName = "LU_E.txt";
		saveToCSV(solve.time, fileName);
	}
}

int main()
{
	int N = 9*100+C*10+D;
	
	Matrix b(N, 1);
	b.initializeBVector(F);	
	
	//cout << "%%% Wektor b %%%" << endl;
	//b.print();

	// ZAD A i B
	zadA_B(N, b);

	// ZAD C
	zadC_D(N, b);

	// ZAD E
	zadE();

	return 0;
}