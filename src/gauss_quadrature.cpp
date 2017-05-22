//
// Created by 季俊哲 on 28/04/2017.
//

#include "gauss_quadrature.h"
#include "basis_function.h"
#include "fe_solution_triangle.h"

/*
 *  gauss quadrature for assemble matrix
 */

double gauss2d_integral_trial_test(vector<vector<double> > vertices,
                                     string basis_type_trial, string basis_type_test,
                                     int basis_index_trial, int basis_index_test,
                                     int trial_derivative_degree_x, int trial_derivative_degree_y,
                                     int test_derivative_degree_x, int test_derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];

    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
//        cout << x << " " << y << endl;
        double trial = basis_function(x, y, vertices, basis_type_trial, basis_index_trial, trial_derivative_degree_x, trial_derivative_degree_y);
        double test = basis_function(x, y, vertices, basis_type_test, basis_index_test, test_derivative_degree_x, test_derivative_degree_y);
//        cout << trial << " " << test << endl;
        quad = quad + gauss_weight[i] * trial * test;
    }
    quad = area * quad;
    return quad;
}

double gauss2d_integral_f_trial_test(double (*func)(double x, double y),
                                     vector<vector<double> > vertices,
                                     string basis_type_trial, string basis_type_test,
                                     int basis_index_trial, int basis_index_test,
                                     int trial_derivative_degree_x, int trial_derivative_degree_y,
                                     int test_derivative_degree_x, int test_derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];

    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * func(x, y)
                                      * basis_function(x, y, vertices, basis_type_trial, basis_index_trial, trial_derivative_degree_x, trial_derivative_degree_y)
                                      * basis_function(x, y, vertices, basis_type_test, basis_index_test, test_derivative_degree_x, test_derivative_degree_y);
    }
    quad = area * quad;
    return quad;
}

double gauss2d_integral_FE_trial_test(VectorXd uh_local, string FE_basis_type, int FE_derivative_degree_x, int FE_derivative_degree_y,
                                     vector<vector<double> > vertices,
                                     string basis_type_trial, int basis_index_trial,
                                     int trial_derivative_degree_x, int trial_derivative_degree_y,
                                     string basis_type_test, int basis_index_test,
                                     int test_derivative_degree_x, int test_derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];

    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * fe_solution_triangle(uh_local, x, y, vertices,FE_basis_type,  FE_derivative_degree_x, FE_derivative_degree_y)
                      * basis_function(x, y, vertices, basis_type_trial, basis_index_trial, trial_derivative_degree_x, trial_derivative_degree_y)
                      * basis_function(x, y, vertices, basis_type_test, basis_index_test, test_derivative_degree_x, test_derivative_degree_y);
    }
    quad = area * quad;
    return quad;
}

double gauss2d_integral_2FE_trial_test(VectorXd uh_local_1, string FE1_basis_type, int FE1_derivative_degree_x, int FE1_derivative_degree_y,
                                       VectorXd uh_local_2, string FE2_basis_type, int FE2_derivative_degree_x, int FE2_derivative_degree_y,
                                       vector <vector<double> > vertices,
                                       string basis_type_trial, int basis_index_trial,
                                       int trial_derivative_degree_x, int trial_derivative_degree_y,
                                       string basis_type_test, int basis_index_test,
                                       int test_derivative_degree_x, int test_derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];

    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * fe_solution_triangle(uh_local_1, x, y, vertices, FE1_basis_type, FE1_derivative_degree_x, FE1_derivative_degree_y)
                      * fe_solution_triangle(uh_local_2, x, y, vertices, FE2_basis_type, FE2_derivative_degree_x, FE2_derivative_degree_y)
                      * basis_function(x, y, vertices, basis_type_trial, basis_index_trial, trial_derivative_degree_x, trial_derivative_degree_y)
                      * basis_function(x, y, vertices, basis_type_test, basis_index_test, test_derivative_degree_x, test_derivative_degree_y);
    }
    quad = area * quad;
    return quad;
}

/*
 *  gauss quadrature for load vector
 */

double gauss2d_integral_test(vector<vector<double>> vertices, int basis_index, string basis_type, int derivative_degree_x, int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];
    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * basis_function(x, y, vertices, basis_type, basis_index, derivative_degree_x, derivative_degree_y);
    }
    quad = quad * area;
    return quad;
}

double gauss2d_integral_f_test(double (*func)(double x, double y),
                               vector<vector<double>> vertices, int basis_index, string basis_type,
                               int derivative_degree_x, int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];
    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * func(x, y)
                                      * basis_function(x, y, vertices, basis_type, basis_index, derivative_degree_x, derivative_degree_y);
    }
    quad = quad * area;
    return quad;
}

double gauss2d_integral_FE_test(VectorXd uh_local, string FE_basis_type, int FE_derivative_degree_x, int FE_derivative_degree_y,
                                vector<vector<double>> vertices,
                                string basis_type, int basis_index,
                                int derivative_degree_x, int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];
    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * fe_solution_triangle(uh_local, x, y, vertices, FE_basis_type, FE_derivative_degree_x, FE_derivative_degree_y)
                      * basis_function(x, y, vertices, basis_type, basis_index, derivative_degree_x, derivative_degree_y);
    }
    quad = quad * area;
    return quad;
}

double gauss2d_integral_2FE_test(VectorXd uh_local_1, string FE1_basis_type, int FE1_derivative_degree_x, int FE1_derivative_degree_y,
                                 VectorXd uh_local_2, string FE2_basis_type, int FE2_derivative_degree_x, int FE2_derivative_degree_y,
                                 vector<vector<double>> vertices,
                                 string basis_type, int basis_index,
                                 int derivative_degree_x, int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];
    double area = fabs(x1 * (y2-y3) + x2 * (y3-y1) + x3 * (y1-y2)) / 2.0;
    vector <vector<double> > gauss_nodes = {{0.33333333333333, 0.33333333333333}, {0.47014206410511, 0.47014206410511},
                                            {0.47014206410511, 0.05971587178977}, {0.05971587178977, 0.47014206410511},
                                            {0.10128650732346, 0.10128650732346}, {0.10128650732346, 0.79742698535309},
                                            {0.79742698535309, 0.10128650732346}};
    vector<double> gauss_weight = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    double quad = 0.0;
    int Ng = (int) gauss_weight.size();
    for (int i = 0; i < Ng; i++) {
        double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
        double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];
        quad = quad + gauss_weight[i] * fe_solution_triangle(uh_local_1, x, y, vertices, FE1_basis_type, FE1_derivative_degree_x, FE1_derivative_degree_y)
                      * fe_solution_triangle(uh_local_2, x, y, vertices, FE2_basis_type, FE2_derivative_degree_x, FE2_derivative_degree_y)
                      * basis_function(x, y, vertices, basis_type, basis_index, derivative_degree_x, derivative_degree_y);
    }
    quad = quad * area;
    return quad;
}