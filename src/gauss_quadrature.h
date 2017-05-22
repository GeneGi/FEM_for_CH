//
// Created by 季俊哲 on 28/04/2017.
//

#ifndef NAVIERSTOKES_EQ_GAUSS_QUADRATURE_H
#define NAVIERSTOKES_EQ_GAUSS_QUADRATURE_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "../lib/Eigen/Core"

using namespace std;
using namespace Eigen;

double gauss2d_integral_trial_test(vector<vector<double> > vertices,
                                     string basis_type_trial, string basis_type_test,
                                     int basis_index_trial, int basis_index_test,
                                     int trial_derivative_degree_x, int trial_derivative_degree_y,
                                     int test_derivative_degree_x, int test_derivative_degree_y);

double gauss2d_integral_f_trial_test(double (*func)(double x, double y),
                                     vector<vector<double> > vertices,
                                     string basis_type_trial, string basis_type_test,
                                     int basis_index_trial, int basis_index_test,
                                     int trial_derivative_degree_x, int trial_derivative_degree_y,
                                     int test_derivative_degree_x, int test_derivative_degree_y);

double gauss2d_integral_FE_trial_test(VectorXd uh_local, string FE_basis_type, int FE_derivative_degree_x, int FE_derivative_degree_y,
                                      vector<vector<double> > vertices,
                                      string basis_type_trial, int basis_index_trial,
                                      int trial_derivative_degree_x, int trial_derivative_degree_y,
                                      string basis_type_test, int basis_index_test,
                                      int test_derivative_degree_x, int test_derivative_degree_y);

double gauss2d_integral_2FE_trial_test(VectorXd uh_local_1, string FE1_basis_type, int FE1_derivative_degree_x, int FE1_derivative_degree_y,
                                       VectorXd uh_local_2, string FE2_basis_type, int FE2_derivative_degree_x, int FE2_derivative_degree_y,
                                       vector <vector<double> > vertices,
                                       string basis_type_trial, int basis_index_trial,
                                       int trial_derivative_degree_x, int trial_derivative_degree_y,
                                       string basis_type_test, int basis_index_test,
                                       int test_derivative_degree_x, int test_derivative_degree_y);

double gauss2d_integral_test(vector<vector<double>> vertices, int basis_index, string basis_type,
                             int derivative_degree_x, int derivative_degree_y);

double gauss2d_integral_f_test(double (*func)(double x, double y),
                               vector<vector<double>> vertices, int basis_index, string basis_type,
                               int derivative_degree_x, int derivative_degree_y);

double gauss2d_integral_FE_test(VectorXd uh_local, string FE_basis_type, int FE_derivative_degree_x, int FE_derivative_degree_y,
                                vector<vector<double>> vertices,
                                string basis_type, int basis_index,
                                int derivative_degree_x, int derivative_degree_y);

double gauss2d_integral_2FE_test(VectorXd uh_local_1, string FE1_basis_type, int FE1_derivative_degree_x, int FE1_derivative_degree_y,
                                 VectorXd uh_local_2, string FE2_basis_type, int FE2_derivative_degree_x, int FE2_derivative_degree_y,
                                 vector<vector<double>> vertices,
                                 string basis_type, int basis_index,
                                 int derivative_degree_x, int derivative_degree_y);

#endif //NAVIERSTOKES_EQ_GAUSS_QUADRATURE_H
