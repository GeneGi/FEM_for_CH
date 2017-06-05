//
// Created by 季俊哲 on 11/05/2017.
//

#include "initial_function.h"
#include <ctime>
#include <cstdlib>



double initial_function_phi(double x, double y) {
    double result;
//    if (x >= -0.5 && x <= 0.5 && y >= -0.5 && y <= 0.5) {
//        result = 1;
//    } else {
//        result = -1;
//    }
    result = -0.05 + ((float)rand() / (float)RAND_MAX) * 0.1;
    return result;
}

double initial_function_u1(double x, double y) {
    return 0;
}

double initial_function_u2(double x, double y) {
    return 0;
}