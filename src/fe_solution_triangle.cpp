//
// Created by 季俊哲 on 27/04/2017.
//

#include <vector>
#include "../lib/Eigen/Sparse"
#include "basis_function.h"

using namespace Eigen;
using namespace std;

double fe_solution_triangle(VectorXd uh_local, double x, double y, vector<vector<double> > vertices,
                            string basis_type, int derivative_degree_x, int derivative_degree_y) {
    int Nlb = uh_local.size();
    double uh_next = 0;
    for (int i = 0; i < Nlb; i++) {
        uh_next += uh_local[i] * basis_function(x, y, vertices, basis_type, i, derivative_degree_x, derivative_degree_y);
    }
    return uh_next;
}
