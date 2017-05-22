//
// Created by 季俊哲 on 10/05/2017.
//

#ifndef CODE_TREAT_BOUNDARY_H
#define CODE_TREAT_BOUNDARY_H

#include "../lib/Eigen/Sparse"
#include <vector>

using namespace std;
using namespace Eigen;

SparseMatrix<double> treat_boundary_A(vector<vector<int> > boundary_nodes, SparseMatrix<double> A);

SparseMatrix<double> treat_boundary_A_three_matrix(
        vector<vector<int> > boundary_nodes_1, int num_of_FE_nodes_1,
        vector<vector<int> > boundary_nodes_2, int num_of_FE_nodes_2,
        vector<vector<int> > boundary_nodes_3,
        SparseMatrix<double> A);

VectorXd treat_boundary_b(double (*boundary_func) (vector<double> omega, double x, double y), vector<double> omega, vector<vector<int> > boundary_nodes, VectorXd b,
                          vector<vector<double> > Pb, double t);

VectorXd treat_boundary_b_three(vector<double> omega,  VectorXd b,
                                double (*boundary_func1) (vector<double> omega, double x, double y), vector<vector<int> > boundary_nodes1, vector<vector<double> > Pb_1, int num_of_FE_nodes_1,
                                double (*boundary_func2) (vector<double> omega, double x, double y), vector<vector<int> > boundary_nodes2, vector<vector<double> > Pb_2, int num_of_FE_nodes_2,
                                double (*boundary_func3) (vector<double> omega, double x, double y), vector<vector<int> > boundary_nodes3, vector<vector<double> > Pb_3);

#endif //CODE_TREAT_BOUNDARY_H
