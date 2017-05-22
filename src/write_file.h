//
// Created by 季俊哲 on 21/05/2017.
//

#ifndef CODE_WRITE_FILE_H
#define CODE_WRITE_FILE_H

#include <fstream>
#include <vector>
#include "../lib/Eigen/Sparse"

using namespace Eigen;
using namespace std;

void sparse2csv(SparseMatrix<double> matrix, string file_name);
void vectro2csv(VectorXd b, string file_name);
void result2dat(VectorXd result, string file_name, vector<vector<double>> Pb, vector<vector<int>> Tb);

#endif //CODE_WRITE_FILE_H
