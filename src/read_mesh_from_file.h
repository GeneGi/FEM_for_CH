//
// Created by 季俊哲 on 24/05/2017.
//

#ifndef CODE_READ_MESH_FROM_FILE_H
#define CODE_READ_MESH_FROM_FILE_H

#include <vector>
#include <fstream>
using namespace std;

vector<vector<double>> read_P_from_file(string file_path);
vector<vector<int>> read_T_from_file(string file_path);
vector<vector<int>> read_boundary_nodes_from_file(string file_path);

#endif //CODE_READ_MESH_FROM_FILE_H
