//
// Created by 季俊哲 on 21/05/2017.
//

#include "write_file.h"

void sparse2csv(SparseMatrix<double> matrix, string file_name) {
    ofstream file(file_name);
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
            if (j == matrix.cols()-1) {
                if (matrix.coeffRef(i, j) == 0) {
                    file << 0;
                } else {
                    file <<  matrix.coeffRef(i, j);
                }
            } else {
                if (matrix.coeffRef(i, j) == 0) {
                    file << 0 << ",";
                } else {
                    file <<  matrix.coeffRef(i, j) << ",";
                }
            }
        }
        file << endl;
    }
    file.close();
}

void vectro2csv(VectorXd b, string file_name) {
    ofstream file(file_name);
    for (int i = 0; i < b.size(); i++) {
        file << b[i] << endl;
    }
    file.close();
}

void result2dat(VectorXd result, string file_name, vector<vector<double>> Pb, vector<vector<int>> Tb)  {
    ofstream file(file_name);
    file << "Variables = X Y phi" << endl;
    file << "Zone N=" << Pb.size() <<  " E=" << Tb.size() << " F=FEPOINT ET=TRIANGLE" << endl;
    for (int i = 0; i < result.size(); i++) {
        file << Pb[i][0] << " " << Pb[i][1] << " " << result[i] << endl;
    }
    for (int i = 0; i < Tb.size(); i++) {
        file << Tb[i][0]+1 << " " << Tb[i][1]+1 << " " << Tb[i][2]+1 << endl;
    }
    file.close();
}
