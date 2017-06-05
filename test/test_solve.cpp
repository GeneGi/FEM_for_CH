//
// Created by 季俊哲 on 18/05/2017.
//
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include "../lib/Eigen/Sparse"
#include "../lib/Eigen/UmfPackSupport"
#include "../lib/Eigen/SPQRSupport"

using namespace std;
using namespace Eigen;

int main()
{
    auto begin = std::chrono::steady_clock::now();
    clock_t time_start = clock();
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
    for (int i = 0; i < 100; i++) {
        ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
//    LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
//    BiCGSTAB<SparseMatrix<double>> solver;  +inf
//    SimplicialLLT<SparseMatrix<double>>  solver; 0
//    SparseLU<SparseMatrix<double>> solver; random
//    SPQR<SparseMatrix<double>> solver;
//    UmfPackLU<SparseMatrix<double>> solver;
//    SPQR<SparseMatrix<double>> solver;
//    CholmodSupernodalLLT<SparseMatrix<double>> solver;
        solver.compute(A);
        VectorXd X = solver.solve(b);
    }
    clock_t time_end = clock();
    cout << "solve time: " << (double)((time_end - time_start)/(CLOCKS_PER_SEC)) << "s" << endl;
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin);
    std::cout << "Cost of method0() is " << elapsed .count() << " milliseconds" << std::endl;
//    for (int m = 0; m < 162; m++) {
//        cout << X[m] << endl;
//    }
//    vectro2csv(X, "result_X.csv");
    return 0;
}
