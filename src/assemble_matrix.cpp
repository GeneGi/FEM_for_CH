//
// Created by 季俊哲 on 28/04/2017.
//

#include "assemble_matrix.h"
#include "gauss_quadrature.h"

#include <iostream>
#include <omp.h>
#include <ctime>
#include <chrono>
#include <ratio>
using namespace std;

SparseMatrix<double> assemble_matrix_A(  vector<vector<double>> P, vector<vector<int>> T,
                                         vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test,
                                         int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                         string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                         string test_basis_type, int test_derivative_degree_x, int test_derivative_degree_y) {
//    MatrixXd StiffnessMatrix(Nb_test, Nb_trial);
    vector<vector<double>> vertices(3, vector<double>(2));
    auto begin = std::chrono::steady_clock::now();
    vector< Triplet<double> > coefficients;
//    clock_t time_start = clock();
    SparseMatrix<double> StiffnessMatrix(Nb_test, Nb_trial);
//    omp_set_num_threads(16);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int alpha = 0; alpha < Nlb_trial; alpha++) {
            for (int beta = 0; beta < Nlb_test; beta++) {
                double temp = gauss2d_integral_trial_test(vertices, trial_basis_type, test_basis_type, alpha, beta,
                                                            trial_derivative_degree_x, trial_derivative_degree_y, test_derivative_degree_x, test_derivative_degree_y);
//                StiffnessMatrix.coeffRef(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
//                StiffnessMatrix(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
                coefficients.push_back(Triplet<double>(Tb_test[i][beta], Tb_trial[i][alpha], temp));
            }
        }
    }
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin);
    std::cout << "Cost of method is " << elapsed .count() << " milliseconds" << std::endl;
//    SparseMatrix<double> result = StiffnessMatrix.sparseView();
    StiffnessMatrix.setFromTriplets(coefficients.begin(), coefficients.end());
//    clock_t time_end = clock();
//    cout << "assemble time: " << (double)((time_end - time_start)/(CLOCKS_PER_SEC)) << "s" << endl;
    return StiffnessMatrix;
//    return result;
}

SparseMatrix<double> assemble_matrix_A_f(double (*func) (double x, double y),
                                       vector<vector<double>> P, vector<vector<int>> T,
                                       vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test,
                                       int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                       string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                       string test_basis_type, int test_derivative_degree_x, int test_derivative_degree_y) {
    SparseMatrix<double> StiffnessMatrix(Nb_test, Nb_trial);
    vector<vector<double>> vertices(3, vector<double>(2));
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int alpha = 0; alpha < Nlb_trial; alpha++) {
            for (int beta = 0; beta < Nlb_test; beta++) {
                double temp = gauss2d_integral_f_trial_test(func, vertices, trial_basis_type, test_basis_type, alpha, beta,
                                                            trial_derivative_degree_x, trial_derivative_degree_y, test_derivative_degree_x, test_derivative_degree_y);
                StiffnessMatrix.coeffRef(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
            }
        }
    }
    return StiffnessMatrix;
}

SparseMatrix<double> assemble_matrix_A_FE(VectorXd un,  string un_basis_type, int un_derivative_degree_x, int un_derivative_degree_y,
                                          vector<vector<double>> P, vector<vector<int>> T,
                                          vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test, vector<vector<int>> Tb_FE,
                                          int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                          string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                          string test_basis_type,int test_derivative_degree_x, int test_derivative_degree_y) {
    SparseMatrix<double> StiffnessMatrix(Nb_test, Nb_trial);
    vector< Triplet<double> > coefficients;
//    MatrixXd StiffnessMatrix(Nb_test, Nb_trial);
    int Nlb = 0;
    if (un_basis_type == "linear") {
        Nlb = 3;
    } else if (un_basis_type == "quadratic") {
        Nlb = 6;
    }
    auto begin = std::chrono::steady_clock::now();
    vector<vector<double>> vertices(3, vector<double>(2));
    VectorXd un_local(Nlb);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int k = 0; k < Nlb; k++) {
            un_local[k] = un[Tb_FE[i][k]];
        }
        for (int alpha = 0; alpha < Nlb_trial; alpha++) {
            for (int beta = 0; beta < Nlb_test; beta++) {
                double temp = gauss2d_integral_FE_trial_test(un_local, un_basis_type, un_derivative_degree_x, un_derivative_degree_y, vertices,
                                                             trial_basis_type, alpha, trial_derivative_degree_x, trial_derivative_degree_y,
                                                             test_basis_type, beta, test_derivative_degree_x, test_derivative_degree_y);
//                StiffnessMatrix.coeffRef(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
//                StiffnessMatrix(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
                coefficients.push_back(Triplet<double>(Tb_test[i][beta], Tb_trial[i][alpha], temp));
            }
        }
    }
//    SparseMatrix<double> result = StiffnessMatrix.sparseView();
    StiffnessMatrix.setFromTriplets(coefficients.begin(), coefficients.end());
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin);
    std::cout << "Cost of method is " << elapsed .count() << " milliseconds" << std::endl;
//    return result;
    return StiffnessMatrix;
}

SparseMatrix<double> assemble_matrix_A_2FE(VectorXd un_1,  string un1_basis_type, int un1_derivative_degree_x, int un1_derivative_degree_y,
                                           VectorXd un_2,  string un2_basis_type, int un2_derivative_degree_x, int un2_derivative_degree_y,
                                          vector<vector<double>> P, vector<vector<int>> T,
                                          vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test, vector<vector<int>> Tb_FE_1, vector<vector<int>> Tb_FE_2,
                                          int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                          string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                          string test_basis_type,int test_derivative_degree_x, int test_derivative_degree_y) {
    SparseMatrix<double> StiffnessMatrix(Nb_test, Nb_trial);
    vector< Triplet<double> > coefficients;
//    MatrixXd StiffnessMatrix(Nb_test, Nb_trial);
    int Nlb_1 = 0;
    int Nlb_2 = 0;
    if (un1_basis_type == "linear") {
        Nlb_1 = 3;
    } else if (un1_basis_type == "quadratic") {
        Nlb_1 = 6;
    }
    if (un2_basis_type == "linear") {
        Nlb_2 = 3;
    } else if (un2_basis_type == "quadratic") {
        Nlb_2 = 6;
    }
    auto begin = std::chrono::steady_clock::now();
    vector<vector<double>> vertices(3, vector<double>(2));
    VectorXd un1_local(Nlb_1), un2_local(Nlb_2);
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int k = 0; k < Nlb_1; k++) {
            un1_local[k] = un_1[Tb_FE_1[i][k]];
        }
        for (int k = 0; k < Nlb_2; k++) {
            un2_local[k] = un_2[Tb_FE_2[k][k]];
        }
        for (int alpha = 0; alpha < Nlb_trial; alpha++) {
            for (int beta = 0; beta < Nlb_test; beta++) {
                double temp = gauss2d_integral_2FE_trial_test(un1_local, un1_basis_type, un1_derivative_degree_x, un1_derivative_degree_y,
                                                              un2_local, un2_basis_type, un2_derivative_degree_x, un2_derivative_degree_y, vertices,
                                                             trial_basis_type, alpha, trial_derivative_degree_x, trial_derivative_degree_y,
                                                             test_basis_type, beta, test_derivative_degree_x, test_derivative_degree_y);
//                StiffnessMatrix.coeffRef(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
//                StiffnessMatrix(Tb_test[i][beta], Tb_trial[i][alpha]) += temp;
                coefficients.push_back(Triplet<double>(Tb_test[i][beta], Tb_trial[i][alpha], temp));
            }
        }
    }
//    SparseMatrix<double> result = StiffnessMatrix.sparseView();
    StiffnessMatrix.setFromTriplets(coefficients.begin(), coefficients.end());
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin);
    std::cout << "Cost of method is " << elapsed .count() << " milliseconds" << std::endl;
//    return result;
    return StiffnessMatrix;
}