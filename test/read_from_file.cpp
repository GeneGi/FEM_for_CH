//
// Created by 季俊哲 on 20/05/2017.
//


#include <iostream>
#include <fstream>
#include <string>

#include "../lib/Eigen/Sparse"

#include "../src/read_mesh_from_file.h"

using namespace Eigen;
using namespace std;

int main()
{
    vector<vector<double>> P = read_P_from_file("../mesh/mesh_P.dat");
    cout << P[0][0] << " " << P[0][1] << endl;
    return 0;
}