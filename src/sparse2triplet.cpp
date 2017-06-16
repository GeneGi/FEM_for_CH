//
// Created by 季俊哲 on 16/06/2017.
//

#include "sparse2triplet.h"


vector<Triplet<double>> sparse2triplet(SparseMatrix<double> sparseMatrix, int row_offset, int col_offset) {
    vector< Triplet<double> > coefficients;
    for (int m = 0; m < sparseMatrix.outerSize(); m++) {
        for (SparseMatrix<double>::InnerIterator it(sparseMatrix, m); it; ++it) {
            coefficients.push_back(Triplet<double>(it.row() + row_offset, it.col() + col_offset, it.value()));
        }
    }
    return coefficients;
}