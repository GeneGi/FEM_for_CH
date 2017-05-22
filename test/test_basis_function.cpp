//
// Created by 季俊哲 on 20/05/2017.
//

#include <iostream>
#include <vector>
#include <cmath>

#include "../src/basis_function.h"

using namespace std;

int main()
{
    vector< vector<double> > vertices = {{-1, -1}, {-1, -0.75}, {-0.75, -1}};
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
    int i = 1;
    double x = x1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + x2 * gauss_nodes[i][0] + x3 * gauss_nodes[i][1];
    double y = y1 * (1 - gauss_nodes[i][0] - gauss_nodes[i][1]) + y2 * gauss_nodes[i][0] + y3 * gauss_nodes[i][1];

    for (int i = 0; i < 6; i++) {
        double result = basis_function(x, y, vertices, "quadratic", i, 0, 0);
        cout << result << endl;
    }

    return 0;
}