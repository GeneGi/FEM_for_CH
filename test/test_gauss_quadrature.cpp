//
// Created by 季俊哲 on 20/05/2017.
//

#include <iostream>
#include <vector>

#include "../src/gauss_quadrature.h"

using namespace std;

int main()
{
    vector< vector<double> > vertices = {{-1, -1}, {-1, -0.75}, {-0.75, -1}};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double result =  gauss2d_integral_trial_test(vertices, "quadratic", "quadratic", i, j, 0, 0, 0, 0);
//            cout << "trial index: " << i << " test index: " << j << " result: " << result << endl;
            cout << result << " ";
        }
        cout << endl;
    }
    return 0;
}