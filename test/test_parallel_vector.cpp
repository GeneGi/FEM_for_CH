//
// Created by 季俊哲 on 09/06/2017.
//

#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;

int main()
{
    vector<vector<double>> v;
#pragma omp parallel for
    for (int i = 0; i < 100000; i++) {
        vector<double> temp = {0, 1, 2, 3, 4};
#pragma omp critical
        v.push_back(temp);
    }
    return 0;

}