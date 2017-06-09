//
// Created by 季俊哲 on 05/05/2017.
//

#include "basis_function.h"

double triangular_reference_basis(double x, double y, string basis_type, int basis_index, int derivative_degree_x, int derivative_degree_y) {
    double result;
    if (basis_type == "linear") {
        if (derivative_degree_x == 0 && derivative_degree_y == 0) {
            if (basis_index == 0) {
                result = 1 - x - y;
            } else if (basis_index == 1) {
                result = x;
            } else if (basis_index == 2) {
                result = y;
            }
        } else if (derivative_degree_x == 1 && derivative_degree_y == 0) {
            if (basis_index == 0) {
                result = -1;
            } else if (basis_index == 1) {
                result = 1;
            } else if (basis_index == 2) {
                result = 0;
            }
        } else if (derivative_degree_x == 0 && derivative_degree_y == 1) {
            if (basis_index == 0) {
                result = -1;
            } else if (basis_index == 1) {
                result = 0;
            } else if (basis_index == 2) {
                result = 1;
            }
        }
    } else if (basis_type == "quadratic") {
        if (derivative_degree_x == 0 && derivative_degree_y == 0) {
            if (basis_index == 0) {
                result = 1 - 3*x - 3*y + 2*x*x + 2*y*y + 4*x*y;
            } else if (basis_index == 1) {
                result = 2*x*x - x;
            } else if (basis_index == 2) {
                result = 2*y*y - y;
            } else if (basis_index == 3) {
                result = 4*x - 4*x*x - 4*x*y;
            } else if (basis_index == 4) {
                result = 4*x*y;
            } else if (basis_index == 5) {
                result = 4*y - 4*y*y - 4*x*y;
            }
        } else if (derivative_degree_x == 1 && derivative_degree_y == 0) {
            if (basis_index == 0) {
                result = -3 + 4*x + 4*y;
            } else if (basis_index == 1) {
                result = 4*x - 1;
            } else if (basis_index == 2) {
                result = 0;
            } else if (basis_index == 3) {
                result = 4 - 8*x - 4*y;
            } else if (basis_index == 4) {
                result = 4*y;
            } else if (basis_index == 5) {
                result = -4*y;
            }
        } else if (derivative_degree_x == 0 && derivative_degree_y == 1) {
            if (basis_index == 0) {
                result = -3 + 4*y + 4*x;
            } else if (basis_index == 1) {
                result = 0;
            } else if (basis_index == 2) {
                result = 4*y - 1;
            } else if (basis_index == 3) {
                result = -4 * x;
            } else if (basis_index == 4) {
                result = 4*x;
            } else if (basis_index == 5) {
                result = 4 - 8*y - 4*x;
            }
        } else if (derivative_degree_x == 2 && derivative_degree_y == 0) {
            if (basis_index == 0) {
                result = 4;
            } else if (basis_index == 1) {
                result = 4;
            } else if (basis_index == 2) {
                result = 0;
            } else if (basis_index == 3) {
                result = -8;
            } else if (basis_index == 4) {
                result = 0;
            } else if (basis_index == 5) {
                result = 0;
            }
        } else if (derivative_degree_x == 0 && derivative_degree_y == 2) {
            if (basis_index == 0) {
                result = 4;
            } else if (basis_index == 1) {
                result = 0;
            } else if (basis_index == 2) {
                result = 4;
            } else if (basis_index == 3) {
                result = 0;
            } else if (basis_index == 4) {
                result = 0;
            } else if (basis_index == 5) {
                result = -8;
            }
        } else if (derivative_degree_x == 1 && derivative_degree_y == 1) {
            if (basis_index == 0) {
                result = 4;
            } else if (basis_index == 1) {
                result = 0;
            } else if (basis_index == 2) {
                result = 0;
            } else if (basis_index == 3) {
                result = -4;
            } else if (basis_index == 4) {
                result = 4;
            } else if (basis_index == 5) {
                result = -4;
            }
        }
    }
    return  result;
}

double basis_function(double x, double y, vector <vector<double> > vertices, string basis_type, int basis_index,
                      int derivative_degree_x, int derivative_degree_y) {
    double x1 = vertices[0][0];
    double y1 = vertices[0][1];
    double x2 = vertices[1][0];
    double y2 = vertices[1][1];
    double x3 = vertices[2][0];
    double y3 = vertices[2][1];

    double J11 = x2 - x1;
    double J12 = x3 - x1;
    double J21 = y2 - y1;
    double J22 = y3 - y1;
    double J = J11 * J22 - J12 * J21;

    double x_hat = (J22 * (x - x1) - J12 * (y - y1)) / J;
    double y_hat = (-J21 * (x - x1) + J11 * (y - y1)) / J;


    double result = 0;
    if (derivative_degree_x == 0 && derivative_degree_y == 0) {
        result = triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 0, 0);
    } else if (derivative_degree_x == 1 && derivative_degree_y == 0) {
        result = (J22 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 1, 0) +
                (-J21) * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 0, 1)) / J;
    } else if (derivative_degree_x == 0 && derivative_degree_y == 1) {
        result = (-J12 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 1, 0) +
                   J11 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 0, 1)) / J;
    } else if (derivative_degree_x == 2 && derivative_degree_y == 0) {
        result = (J22 * J22 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 2, 0) +
                  J21 * J21 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 0, 2) +
                (-2*J21 * J22) * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 1, 1)) / (J * J);
    } else if (derivative_degree_x == 0 && derivative_degree_y == 2) {
        result = (J12 * J12 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 2, 0) +
                  J11 * J11 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 0, 2) +
                  (-2*J11 * J12) * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 1, 1)) / (J * J);
    } else if (derivative_degree_x == 1 && derivative_degree_y == 1) {
        result = (-J22 * J12 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 2, 0) +
                  -J21 * J11 * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 0, 2) +
                  (J21 * J12 + J11 * J22) * triangular_reference_basis(x_hat, y_hat, basis_type, basis_index, 1, 1)) / (J * J);
    }

    return result;
}