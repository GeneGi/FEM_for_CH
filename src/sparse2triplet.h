//
// Created by 季俊哲 on 16/06/2017.
//

#ifndef CODE_SPARSE2TRIPLET_H
#define CODE_SPARSE2TRIPLET_H

#include <vector>
#include "../lib/Eigen/Sparse"

using namespace std;
using namespace Eigen;

vector<Triplet<double>> sparse2triplet(SparseMatrix<double> sparseMatrix, int row_offset, int col_offset);

#endif //CODE_SPARSE2TRIPLET_H
