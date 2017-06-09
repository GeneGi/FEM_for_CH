//
// Created by 季俊哲 on 09/06/2017.
//

#include "fix_pressure.h"
#include "boundary_function.h"

SparseMatrix<double> fix_pressure_A (SparseMatrix<double> A, int num_nodes_before, int fix_p_index) {
    for (int i = 0; i < A.cols(); i++) {
        A.coeffRef(num_nodes_before + fix_p_index, i) = 0;
    }
    A.coeffRef(num_nodes_before + fix_p_index, num_nodes_before + fix_p_index) = 1;
    return A;
}

VectorXd fix_pressure_b(double (*boundary_func) (vector<double> omega, double x, double y), VectorXd b, int num_nodes_before, int fix_p_index, vector<double> omega, vector<vector<double>> M_basis_p) {
    b[num_nodes_before + fix_p_index] = boundary_func(omega, M_basis_p[fix_p_index][0], M_basis_p[fix_p_index][1]);
    return b;
}