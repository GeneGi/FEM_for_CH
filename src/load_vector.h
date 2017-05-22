//
// Created by 季俊哲 on 05/05/2017.
//

#ifndef NAVIERSTOKES_EQ_LOAD_VECTOR_H
#define NAVIERSTOKES_EQ_LOAD_VECTOR_H

#include <vector>
#include "../lib/Eigen/Core"

using namespace std;
using namespace Eigen;

VectorXd load_vector_b(vector<vector<double>> P, vector<vector<double>> Pb, vector<vector<int>> T, vector<vector<int>> Tb,
                       int N, int Nb, int Nlb, string basis_type, int derivative_degree_x, int derivative_degree_y);

VectorXd load_vector_f_b(double (*func)(double x, double y),
                         vector<vector<double>> P, vector<vector<int>> T, vector<vector<int>> Tb,
                         int Nlb, int N, int Nb, string basis_type, int derivative_degree_x, int derivative_degree_y);

VectorXd load_vector_FE_b(VectorXd un,  string un_basis_type, int un_derivative_degree_x, int un_derivative_degree_y,
                          vector<vector<double>> P, vector<vector<int>> T, vector<vector<int>> Tb, vector<vector<int>> Tb_FE,
                          int N, int Nb, int Nlb, string basis_type, int derivative_degree_x, int derivative_degree_y);

VectorXd load_vector_2FE_b(VectorXd un1,  string un1_basis_type, int un1_derivative_degree_x, int un1_derivative_degree_y,
                           VectorXd un2,  string un2_basis_type, int un2_derivative_degree_x, int un2_derivative_degree_y,
                           vector<vector<double>> P, vector<vector<int>> T, vector<vector<int>> Tb, vector<vector<int>> Tb_FE_1, vector<vector<int>> Tb_FE_2,
                           int N, int Nb_test, int Nlb_test, string basis_type, int derivative_degree_x, int derivative_degree_y);

#endif //NAVIERSTOKES_EQ_LOAD_VECTOR_H
