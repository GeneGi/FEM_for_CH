//
// Created by 季俊哲 on 18/05/2017.
//
#include <iostream>
#include <string>
#include "../lib/Eigen/Sparse"

using namespace std;
using namespace Eigen;

int main()
{
    SparseMatrix<double> A(2, 3);

    A.coeffRef(0, 0) = 1;
    A.coeffRef(0, 1) = 1;
    A.coeffRef(0, 2) = 1;
    A.coeffRef(1, 0) = 1;
    A.coeffRef(1, 1) = 1;
    A.coeffRef(1, 2) = 1;

    VectorXd b(2);
    b << 3, 3;

    LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
    solver.compute(A);
    VectorXd X = solver.solve(b);

    for (int i = 0; i < 3; i++) {
        cout << X[i] << endl;
    }
    string name;
    name = (string)"a" + (string)"b";
    return 0;
}
