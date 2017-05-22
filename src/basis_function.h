//
// Created by 季俊哲 on 05/05/2017.
//

#ifndef NAVIERSTOKES_EQ_BASIS_FUNCTION_H
#define NAVIERSTOKES_EQ_BASIS_FUNCTION_H

#include <vector>
#include <string>
#include <iostream>
using namespace std;

double triangular_reference_basis(double x, double y, string basis_type, int basis_index, int derivative_degree_x, int derivative_degree_y);
double basis_function(double x, double y, vector <vector<double> > vertices, string basis_type, int basis_index,
                      int derivative_degree_x, int derivative_degree_y);

#endif //NAVIERSTOKES_EQ_BASIS_FUNCTION_H
