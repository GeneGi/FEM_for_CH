//
// Created by 季俊哲 on 10/05/2017.
//

#include "treat_boundary.h"

SparseMatrix<double> treat_boundary_A(vector<vector<int> > boundary_nodes,
                                      SparseMatrix<double> A) {
    long nbn = boundary_nodes.size();
    for (int i = 0; i < nbn; i++) {
        if (boundary_nodes[i][0] == -1) {
            int k = boundary_nodes[i][1];
            for (int j = 0; j< A.innerSize(); j++) {
                A.coeffRef(k, j) = 0;
            }
            A.coeffRef(k, k) = 1;
        }
    }
    return A;
}

SparseMatrix<double> treat_boundary_A_three_matrix(
                                      vector<vector<int> > boundary_nodes_1, int num_of_FE_nodes_1,
                                      vector<vector<int> > boundary_nodes_2, int num_of_FE_nodes_2,
                                      vector<vector<int> > boundary_nodes_3,
                                      SparseMatrix<double> A) {
    long nbn_1 = boundary_nodes_1.size();
    for (int i = 0; i < nbn_1; i++) {
//        if (boundary_nodes_1[i][0] == -1) {
            int k = boundary_nodes_1[i][1];
            for (int j = 0; j< A.cols(); j++) {
                A.coeffRef(k, j) = 0;
            }
            A.coeffRef(k, k) = 1;
//        }
    }

    long nbn_2 = boundary_nodes_2.size();
    for (int i = 0; i < nbn_2; i++) {
//        if (boundary_nodes_2[i][0] == -1) {
            int k = boundary_nodes_2[i][1];
            for (int j = 0; j< A.cols(); j++) {
                A.coeffRef(k+num_of_FE_nodes_1, j) = 0;
            }
            A.coeffRef(k+num_of_FE_nodes_1, k+num_of_FE_nodes_1) = 1;
//        }
    }

    for (int i = 0; i < A.cols(); i++) {
        A.coeffRef(num_of_FE_nodes_1 + num_of_FE_nodes_2, i) = 0;
    }
    A.coeffRef(num_of_FE_nodes_1 + num_of_FE_nodes_2, num_of_FE_nodes_1 + num_of_FE_nodes_2) = 1;
//    long nbn_3 = boundary_nodes_3.size();
//    for (int i = 0; i < nbn_3; i++) {
//        if (boundary_nodes_3[i][0] == -1) {
//            int k = boundary_nodes_3[i][1];
//            for (int j = 0; j< A.cols(); j++) {
//                A.coeffRef(k+num_of_FE_nodes_1+num_of_FE_nodes_2, j) = 0;
//            }
//            A.coeffRef(k+num_of_FE_nodes_1+num_of_FE_nodes_2, k+num_of_FE_nodes_1+num_of_FE_nodes_2) = 1;
//        }
//    }
    return A;
}

VectorXd treat_boundary_b(double (*boundary_func) (vector<double> omega, double x, double y), vector<double> omega, vector<vector<int> > boundary_nodes, VectorXd b,
                              vector<vector<double> > Pb) {
    long nbn = boundary_nodes.size();
    for (int i = 0; i < nbn; i++) {
        if (boundary_nodes[i][0] == -1) {
            int k = boundary_nodes[i][1];
            b[k] = boundary_func(omega, Pb[k][0], Pb[k][1]);
        }
    }
    return b;
}

VectorXd treat_boundary_b_three(vector<double> omega,  VectorXd b,
                                double (*boundary_func1) (vector<double> omega, double x, double y), vector<vector<int> > boundary_nodes1, vector<vector<double> > Pb_1, int num_of_FE_nodes_1,
                                double (*boundary_func2) (vector<double> omega, double x, double y), vector<vector<int> > boundary_nodes2, vector<vector<double> > Pb_2, int num_of_FE_nodes_2,
                                double (*boundary_func3) (vector<double> omega, double x, double y), vector<vector<int> > boundary_nodes3, vector<vector<double> > Pb_3)
{
    long nbn_1 = boundary_nodes1.size();

    for (int i = 0; i < nbn_1; i++) {
//        if (boundary_nodes1[i][0] == -1) {
            int k = boundary_nodes1[i][1];
            b[k] = boundary_func1(omega, Pb_1[k][0], Pb_1[k][1]);
//        }
    }

    long nbn_2 = boundary_nodes2.size();
    for (int i = 0; i < nbn_2; i++) {
//        if (boundary_nodes2[i][0] == -1) {
            int k = boundary_nodes1[i][1];
            b[k+num_of_FE_nodes_1] = boundary_func2(omega, Pb_2[k][0], Pb_2[k][1]);
//        }
    }

    b[num_of_FE_nodes_1 + num_of_FE_nodes_2] = boundary_func3(omega, Pb_3[0][0], Pb_3[0][1]);
//    for (int i = 0; i < nbn_3; i++) {
//        if (boundary_nodes3[i][0] == -1) {
//            int k = boundary_nodes3[i][1];
//            b[k+num_of_FE_nodes_1+num_of_FE_nodes_2] = boundary_func3(omega, Pb_3[k][0], Pb_3[k][1]);
//        }
//    }
    return b;
}

