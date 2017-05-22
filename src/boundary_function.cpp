//
// Created by 季俊哲 on 11/05/2017.
//

#include "boundary_function.h"

double boundary_function_phi(vector<double> omega, double x, double y) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double result = 0;

    if (x == left) {
        result = 0;
    } else if (x == right) {
        result = 0;
    } else if (y == bottom) {
        result = 0;
    } else if (y == top) {
        result = 0;
    }
    return result;
}

double boundary_function_u1(vector<double> omega, double x, double y) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double result = 0;

    if (x == left) {
        result = 0;
    } else if (x == right) {
        result = 0;
    } else if (y == bottom) {
        result = 0;
    } else if (y == top) {
        result = 0;
    }
    return result;
}

double boundary_function_u2(vector<double> omega, double x, double y) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];
    double result = 0;

    if (x == left) {
        result = 0;
    } else if (x == right) {
        result = 0;
    } else if (y == bottom) {
        result = 0;
    } else if (y == top) {
        result = 0;
    }
    return result;
}

double boundary_function_p(vector<double> omega, double x, double y) {
    double left = omega[0];
    double right = omega[1];
    double bottom = omega[2];
    double top = omega[3];

    double result = 0;

    if (x == left) {
        result = 0;
    } else if (x == right) {
        result = 0;
    } else if (y == bottom) {
        result = 0;
    } else if (y == top) {
        result = 0;
    }
    return result;
}