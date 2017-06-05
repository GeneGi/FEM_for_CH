//
// Created by 季俊哲 on 24/05/2017.
//

#include "read_mesh_from_file.h"

vector<vector<double>> read_P_from_file(string file_path) {
    ifstream in(file_path);
    int N;
    in >> N;
    vector<vector<double>> P(N, vector<double>(2));

    for (int i = 0; i < N; i++) {
        in >> P[i][0];
        in >> P[i][1];
    }
    return P;
}

vector<vector<int>> read_T_from_file(string file_path) {
    ifstream in(file_path);
    int N;
    int Nlb;
    in >> N >> Nlb;
    vector<vector<int>> T(N, vector<int>(Nlb));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < Nlb; j++) {
            in >> T[i][j];
        }
    }
    return T;
}

vector<vector<int>> read_boundary_nodes_from_file(string file_path) {
    ifstream in(file_path);
    int N;
    in >> N;
    vector<vector<int>> boundary_nodes(N, vector<int>(2, -1));
    for (int i = 0; i < N; i++) {
        in >> boundary_nodes[i][0] >> boundary_nodes[i][1];
    }
    return boundary_nodes;
}