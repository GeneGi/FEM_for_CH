//
// Created by 季俊哲 on 27/04/2017.
//

#include <iostream>
#include <omp.h>
#include "generate_info_matrix.h"

using namespace std;

vector <vector<double>> generate_P_triangle(vector<double> omega, vector<double> h, string basis_type) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double hx, hy;

    if (basis_type == "linear") {
        hx = h[0];
        hy = h[1];
    } else if (basis_type == "quadratic") {
        hx = h[0] / 2;
        hy = h[1] / 2;
    } else {
        hx = 0;
        hy = 0;
    }
    int Nx = (int) ((right - left) / hx);
    int Ny = (int) ((top - bottom) / hy);
    int Nb = (Nx + 1) * (Ny + 1);

    vector <vector<double>> P(Nb, vector<double>(2));

    for (int i = 0; i < Nx + 1; i++) {
        for (int j = 0; j < Ny + 1; j++) {
            P[i * (Ny + 1) + j][0] = left + i * hx;
            P[i * (Ny + 1) + j][1] = bottom + j * hy;
        }
    }

    return P;
}

vector <vector<int>> generate_T_triangle(vector<double> omega, vector<double> h, string basis_type) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double hx = h[0];
    double hy = h[1];

    int Nx = (int) ((right - left) / hx);
    int Ny = (int) ((top - bottom) / hy);
    int N = 2 * Nx * Ny;

    if (basis_type == "linear") {
        int Nlb = 3;

        vector<vector<int>> T((unsigned long) N, vector<int>(Nlb));
        vector<vector<int>> temp((unsigned long) (Nx + 1), vector<int>((unsigned long) (Ny + 1)));
        for (int i = 0; i < Nx+1; i++) {
            for (int j = 0; j < Ny+1; j++) {
                temp[i][j] = i * (Ny+1) + j;
            }
        }

        int row, column;
        for (int i = 0; i < Nx*Ny; i++) {
            if ((i+1) % Ny == 0) {
                row = Ny - 1;
                column = i / Ny;
            } else {
                row = i % Ny;
                column = i / Ny;
            }
            T[2*i][0] = temp[column][row];
            T[2*i][1] = temp[column+1][row];
            T[2*i][2] = temp[column][row+1];

            T[2*i+1][0] = temp[column][row+1];
            T[2*i+1][1] = temp[column+1][row];
            T[2*i+1][2] = temp[column+1][row+1];
        }
        return T;
    } else if (basis_type == "quadratic") {
        int Nlb = 6;

        vector<vector<int>> T((unsigned long) N, vector<int>(Nlb));
        vector<vector<int>> temp((unsigned long) (2*Nx + 1), vector<int>((unsigned long) (2*Ny + 1)));
        for (int i = 0; i < 2*Nx+1; i++) {
            for (int j = 0; j < 2*Ny+1; j++) {
                temp[i][j] = i * (2*Ny+1) + j;
            }
        }
        int row, column;
        for (int i = 0; i < Nx*Ny; i++) {
            if ((i+1) % Ny == 0) {
                row = Ny - 1;
                column = i / Ny;
            } else {
                row = i % Ny;
                column = i / Ny;
            }
            T[2*i][0] = temp[2*column][2*row];
            T[2*i][1] = temp[2*column + 2][2*row];
            T[2*i][2] = temp[2*column][2*row + 2];
            T[2*i][3] = temp[2*column + 1][2*row];
            T[2*i][4] = temp[2*column + 1][2*row + 1];
            T[2*i][5] = temp[2*column][2*row + 1];

            T[2*i+1][0] = temp[2*column][2*row + 2];
            T[2*i+1][1] = temp[2*column + 2][2*row];
            T[2*i+1][2] = temp[2*column + 2][2*row + 2];
            T[2*i+1][3] = temp[2*column + 1][2*row + 1];
            T[2*i+1][4] = temp[2*column + 2][2*row + 1];
            T[2*i+1][5] = temp[2*column + 1][2*row + 2];
        }
        return T;
    }
}

vector<vector<int>> generate_boundary_nodes(vector<double> omega, vector<double> h, string basis_type) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double hx = h[0];
    double hy = h[1];

    int Nx = (int) ((right - left) / hx);
    int Ny = (int) ((top - bottom) / hy);

    if (basis_type == "quadratic") {
        Nx = 2 * Nx;
        Ny = 2 * Ny;
    }
    int nbn = 2 * (Nx+Ny);
    vector<vector<int>> boundary_nodes(nbn, vector<int>(2, -1));
    // left side
    for (int i = 0; i < Nx; i++) {
        boundary_nodes[i][1] = i * (Ny+1);
    }
    // top side
    for (int i = Nx; i < Nx+Ny; i++) {
        boundary_nodes[i][1] = Nx * Ny + i;
    }
    // right side
    for (int i = Nx+Ny; i < 2*Nx+Ny; i++) {
        boundary_nodes[i][1] = (2*Nx + Ny + 1 - i) * (Ny + 1) - 1;
    }
    // bottom side
    for (int i = 2*Nx+Ny; i < nbn; i++) {
        boundary_nodes[i][1] = 2*Nx + 2*Ny - i;
    }
    return boundary_nodes;
}
