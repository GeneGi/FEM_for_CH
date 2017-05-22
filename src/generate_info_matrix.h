//
// Created by 季俊哲 on 27/04/2017.
//

#ifndef NAVIERSTOKES_EQ_GENERATE_INFO_MATRIX_H
#define NAVIERSTOKES_EQ_GENERATE_INFO_MATRIX_H

#include <vector>
#include <string>
using  namespace std;

vector <vector<double>> generate_P_triangle(vector<double> omega, vector<double> h, string basis_type);
vector <vector<int>> generate_T_triangle(vector<double> omega, vector<double> h, string basis_type);
vector<vector<int>> generate_boundary_nodes(vector<double> omega, vector<double> h, string basis_type);

#endif //NAVIERSTOKES_EQ_GENERATE_INFO_MATRIX_H
