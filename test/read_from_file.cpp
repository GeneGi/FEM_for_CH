//
// Created by 季俊哲 on 20/05/2017.
//


#include <iostream>
#include <fstream>
#include <string>

#include "../lib/Eigen/Sparse"
#include "../lib/Eigen/IterativeLinearSolvers"
#include "../src/write_file.h"

using namespace Eigen;
using namespace std;

int main()
{
    SparseMatrix<double> A(162, 162);
    ifstream file_A("../data/A001.txt");
    double value;
    int i=0, j=0;
    while (file_A >> value) {
        A.coeffRef(i, j) = value;
        j++;
        if (j == 162) {
            i++;
            j=0;
        }
    }
    VectorXd b(162);
    ifstream file_b("../data/b0.01.txt");
    int k = 0;
    while (file_b >> value) {
        b[k] = value;
        k++;
    }
    ConjugateGradient<SparseMatrix<double>, Upper> solver;
//    LeastSquaresConjugateGradient<SparseMatrix<double>> solver; x random
//    BiCGSTAB<SparseMatrix<double>> solver;  +inf
//    SimplicialLLT<SparseMatrix<double>>  solver; 0
//    SparseLU<SparseMatrix<double>> solver; random
    solver.compute(A);
    VectorXd X = solver.solve(b);
    for (int m = 0; m < 162; m++) {
        cout << X[m] << endl;
    }
//    vectro2csv(X, "result_X.csv");

    return 0;
}