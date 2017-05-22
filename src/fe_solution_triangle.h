//
// Created by 季俊哲 on 27/04/2017.
//

#ifndef NS_EQUATION_FE_SOLUTION_TRIANGLE_H
#define NS_EQUATION_FE_SOLUTION_TRIANGLE_H

#include "../lib/Eigen/Sparse"

using namespace Eigen;
using namespace std;

double fe_solution_triangle(VectorXd uh_local, double x, double y, vector<vector<double> > vertices,
                            string basis_type, int derivative_degree_x, int derivative_degree_y);

#endif //NS_EQUATION_FE_SOLUTION_TRIANGLE_H
