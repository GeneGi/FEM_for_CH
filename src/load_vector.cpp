//
// Created by 季俊哲 on 05/05/2017.
//

#include "load_vector.h"
#include "gauss_quadrature.h"

VectorXd load_vector_b(vector<vector<double>> P, vector<vector<double>> Pb, vector<vector<int>> T, vector<vector<int>> Tb,
                       int N, int Nb, int Nlb, string basis_type, int derivative_degree_x, int derivative_degree_y) {
    VectorXd right_side(Nb);
    for (int i=0; i < Nb; i++) {
        right_side[i] = 0;
    }
    vector<vector<double>> vertices(3, vector<double>(2));
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int beta = 0; beta < Nlb; beta++) {
            double temp = gauss2d_integral_test(vertices, beta, basis_type, derivative_degree_x, derivative_degree_y);
            right_side(Tb[i][beta]) += temp;
        }
    }
    return right_side;
}

VectorXd load_vector_f_b(double (*func)(double x, double y),
                       vector<vector<double>> P, vector<vector<int>> T, vector<vector<int>> Tb,
                         int Nlb, int N, int Nb, string basis_type, int derivative_degree_x, int derivative_degree_y) {
    VectorXd right_side(Nb);
    for (int i=0; i < Nb; i++) {
        right_side[i] = 0;
    }
    vector<vector<double>> vertices(3, vector<double>(2));
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int beta = 0; beta < Nlb; beta++) {
            double temp = gauss2d_integral_f_test(func, vertices, beta, basis_type, derivative_degree_x, derivative_degree_y);
            right_side(Tb[i][beta]) += temp;
        }
    }
    return right_side;
}

VectorXd load_vector_FE_b(VectorXd un,  string un_basis_type, int un_derivative_degree_x, int un_derivative_degree_y,
                       vector<vector<double>> P, vector<vector<int>> T, vector<vector<int>> Tb, vector<vector<int>> Tb_FE,
                       int N, int Nb_test, int Nlb_test, string basis_type, int derivative_degree_x, int derivative_degree_y) {
    VectorXd right_side(Nb_test);
    int Nlb = 0;
    for (int i=0; i < Nb_test; i++) {
        right_side[i] = 0;
    }
    if (un_basis_type == "linear") {
        Nlb = 3;
    } else if (un_basis_type == "quadratic") {
        Nlb = 6;
    }
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
        for (int beta = 0; beta < Nlb_test; beta++) {
            double temp = gauss2d_integral_FE_test(un_local, un_basis_type, un_derivative_degree_x, un_derivative_degree_y,
                                                   vertices,  basis_type, beta, derivative_degree_x, derivative_degree_y);
            right_side(Tb[i][beta]) += temp;
        }
    }
    return right_side;
}

VectorXd load_vector_2FE_b(VectorXd un1,  string un1_basis_type, int un1_derivative_degree_x, int un1_derivative_degree_y,
                           VectorXd un2,  string un2_basis_type, int un2_derivative_degree_x, int un2_derivative_degree_y,
                          vector<vector<double>> P, vector<vector<int>> T, vector<vector<int>> Tb, vector<vector<int>> Tb_FE_1, vector<vector<int>> Tb_FE_2,
                          int N, int Nb_test, int Nlb_test, string basis_type, int derivative_degree_x, int derivative_degree_y) {
    VectorXd right_side(Nb_test);
    for (int i=0; i < Nb_test; i++) {
        right_side[i] = 0;
    }
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

    vector<vector<double>> vertices(3, vector<double>(2));
    VectorXd un1_local(Nlb_1), un2_local(Nlb_2);

    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 3; k++) {
            vertices[k][0] = P[T[i][k]][0];
            vertices[k][1] = P[T[i][k]][1];
        }
        for (int k = 0; k < Nlb_1; k++) {
            un1_local[k] = un1[Tb_FE_1[i][k]];
        }
        for (int k = 0; k < Nlb_2; k++) {
            un2_local[k] = un2[Tb_FE_2[i][k]];
        }
        for (int beta = 0; beta < Nlb_test; beta++) {
            double temp = gauss2d_integral_2FE_test(un1_local, un1_basis_type, un1_derivative_degree_x, un1_derivative_degree_y,
                                                   un2_local, un2_basis_type, un2_derivative_degree_x, un2_derivative_degree_y,
                                                   vertices,  basis_type, beta, derivative_degree_x, derivative_degree_y);
            right_side(Tb[i][beta]) += temp;
        }
    }
    return right_side;
}