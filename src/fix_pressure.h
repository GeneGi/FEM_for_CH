//
// Created by 季俊哲 on 09/06/2017.
//

#ifndef CODE_FIX_PRESSURE_H
#define CODE_FIX_PRESSURE_H

#include "../lib/Eigen/Sparse"
#include <vector>

using namespace std;
using namespace Eigen;

SparseMatrix<double> fix_pressure_A (SparseMatrix<double> A, int num_nodes_before, int fix_p_index);
VectorXd fix_pressure_b(double (*boundary_func) (vector<double> omega, double x, double y), VectorXd b, int num_nodes_before,  int fix_p_index,  vector<double> omega, vector<vector<double>> M_basis_p);


#endif //CODE_FIX_PRESSURE_H
