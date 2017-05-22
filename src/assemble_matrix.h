//
// Created by 季俊哲 on 28/04/2017.
//

#ifndef NAVIERSTOKES_EQ_ASSEMBLE_MATRIX_H
#define NAVIERSTOKES_EQ_ASSEMBLE_MATRIX_H

#include "../lib/Eigen/Sparse"
#include <vector>

using namespace std;
using namespace Eigen;

SparseMatrix<double> assemble_matrix_A(  vector<vector<double>> P, vector<vector<int>> T,
                                         vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test,
                                         int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                         string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                         string test_basis_type, int test_derivative_degree_x, int test_derivative_degree_y);

SparseMatrix<double> assemble_matrix_A_f(double (*func) (double x, double y),
                                         vector<vector<double>> P, vector<vector<int>> T,
                                         vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test,
                                         int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                         string trial_basis_type, string test_basis_type,
                                         int trial_derivative_degree_x, int trial_derivative_degree_y,
                                         int test_derivative_degree_x, int test_derivative_degree_y);

SparseMatrix<double> assemble_matrix_A_FE(VectorXd un,  string un_basis_type, int un_derivative_degree_x, int un_derivative_degree_y,
                                          vector<vector<double>> P, vector<vector<int>> T,
                                          vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test, vector<vector<int>> Tb_FE,
                                          int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                          string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                          string test_basis_type,int test_derivative_degree_x, int test_derivative_degree_y);

SparseMatrix<double> assemble_matrix_A_2FE(VectorXd un_1,  string un1_basis_type, int un1_derivative_degree_x, int un1_derivative_degree_y,
                                           VectorXd un_2,  string un2_basis_type, int un2_derivative_degree_x, int un2_derivative_degree_y,
                                           vector<vector<double>> P, vector<vector<int>> T,
                                           vector<vector<int>> Tb_trial, vector<vector<int>> Tb_test, vector<vector<int>> Tb_FE_1, vector<vector<int>> Tb_FE_2,
                                           int N, int Nb_trial, int Nb_test, int Nlb_trial, int Nlb_test,
                                           string trial_basis_type, int trial_derivative_degree_x, int trial_derivative_degree_y,
                                           string test_basis_type,int test_derivative_degree_x, int test_derivative_degree_y);

#endif //NAVIERSTOKES_EQ_ASSEMBLE_MATRIX_H
