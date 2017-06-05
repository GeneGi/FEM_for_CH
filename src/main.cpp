#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctime>
#include <chrono>

#include "../lib/Eigen/Sparse"

#include "generate_info_matrix.h"
#include "assemble_matrix.h"
#include "initial_function.h"
#include "load_vector.h"
#include "treat_boundary.h"
#include "boundary_function.h"
#include "write_file.h"
#include "read_mesh_from_file.h"

using namespace std;
using namespace Eigen;

int main() {
    initParallel();
    srand(time(NULL));
    omp_set_num_threads(8);
    clock_t time_start = clock();
    auto begin = std::chrono::steady_clock::now();
    cout << "input parameters and caculate intermediate variables:" << endl;
    /*
     *   parameters
     */
    vector<double> Omega = {-1.0, 1.0, -1.0, 1.0};
    vector<double> h = {1.0/32, 1.0/32};
//    vector<double> h = {1.0/16, 1.0/16};
//    vector<double> h = {0.25, 0.25};

    double gamma = 0.01;
    double eta = 0.02;
    double rho_0 = 1;
    double mu_0 = 1;
    double lambda = 1;

    string basis_type_phi = "linear";
    string basis_type_omega = "linear";
    string basis_type_u = "quadratic";
    string basis_type_p = "linear";

    double t_min = 0;
    double t_max = 4;
    double dt = 0.001;

    double left = Omega[0];
    double right = Omega[1];
    double bottom = Omega[2];
    double top = Omega[3];
    double hx = h[0];
    double hy = h[1];

    /*
     *  intermediate variables
     */

    int Nx_partition = (int) ((right - left) / hx);
    int Ny_partition = (int) ((top - bottom) / hy);
    int Nt = (int) ((t_max - t_min) / dt);

    int Nx_basis_phi = Nx_partition;
    int Ny_basis_phi = Ny_partition;
    int Nx_basis_omega = Nx_partition;
    int Ny_basis_omega = Ny_partition;
    int Nx_basis_u = Nx_partition * 2;
    int Ny_basis_u = Ny_partition * 2;
    int Nx_basis_p = Nx_partition;
    int Ny_basis_p = Ny_partition;

    /*
     *  generate info matrix and boundary nodes
     */
//    clock_t time_variables = clock();
//    cout << ((time_variables - time_start) / CLOCKS_PER_SEC) << endl;
//    cout << "generate info matrix and boundary nodes:";
    
    vector< vector<double> > M_partition = generate_P_triangle(Omega, h, "linear");
    vector< vector<int> > T_partition = generate_T_triangle(Omega, h, "linear");
    vector< vector<double> > M_basis_phi = M_partition;
    vector< vector<int> > T_basis_phi = T_partition;
    vector< vector<double> > M_basis_omega = M_partition;
    vector< vector<int> > T_basis_omega = T_partition;
    vector< vector<double> > M_basis_u = generate_P_triangle(Omega, h, basis_type_u);
    vector< vector<int> > T_basis_u = generate_T_triangle(Omega, h, basis_type_u);
    vector< vector<double> > M_basis_p = M_partition;
    vector< vector<int> > T_basis_p = T_partition;

    vector<vector<int>> boundary_nodes_phi = generate_boundary_nodes(Omega, h, basis_type_phi);
    vector<vector<int>> boundary_nodes_omega = generate_boundary_nodes(Omega, h, basis_type_omega);
    vector<vector<int>> boundary_nodes_u = generate_boundary_nodes(Omega, h, basis_type_u);
    vector<vector<int>> boundary_nodes_p= generate_boundary_nodes(Omega, h, basis_type_p);

    int num_of_elements = 2 * Nx_partition * Ny_partition;
    int num_of_FE_nodes_phi = (Nx_basis_phi + 1) * (Ny_basis_phi + 1);
    int num_of_FE_nodes_omega = (Nx_basis_omega + 1) * (Ny_basis_omega + 1);
    int num_of_FE_nodes_u = (Nx_basis_u + 1) * (Ny_basis_u + 1);
    int num_of_FE_nodes_p = (Nx_basis_p + 1) * (Ny_basis_p + 1);

    /*
     *  read mesh from file
     */

//    vector<vector<double>> M_partition = read_P_from_file("../mesh/mesh_P.dat");
//    vector<vector<int>> T_partition = read_T_from_file("../mesh/mesh_T.dat");
//    vector<vector<double>> M_basis_phi = M_partition;
//    vector<vector<int>> T_basis_phi = T_partition;
//    vector<vector<double>> M_basis_omega = M_partition;
//    vector<vector<int>> T_basis_omega = T_partition;
//    vector<vector<double>> M_basis_u = read_P_from_file("../mesh/FE_P.dat");
//    vector<vector<int>> T_basis_u = read_T_from_file("../mesh/FE_T.dat");
//    vector<vector<double>> M_basis_p = M_partition;
//    vector<vector<int>> T_basis_p = T_partition;
//
//    vector<vector<int>> boundary_nodes_phi = read_boundary_nodes_from_file("../mesh/mesh_boundary_nodes.dat");
//    vector<vector<int>> boundary_nodes_omega = boundary_nodes_phi;
//    vector<vector<int>> boundary_nodes_u = read_boundary_nodes_from_file("../mesh/FE_boundary_nodes.dat");
//    vector<vector<int>> boundary_nodes_p = boundary_nodes_phi;
//
//    int num_of_elements = T_partition.size();
//    int num_of_FE_nodes_phi = M_basis_phi.size();
//    int num_of_FE_nodes_omega = M_basis_omega.size();
//    int num_of_FE_nodes_u = M_basis_u.size();
//    int num_of_FE_nodes_p = M_basis_p.size();

    cout << "initial values:" << endl;
    


    int num_of_local_basis_phi = 3;
    int num_of_local_basis_omega = 3;
    int num_of_local_basis_u = 6;
    int num_of_local_basis_p = 3;

    // create result object
    vector<VectorXd> phi(Nt+1, VectorXd(num_of_FE_nodes_phi));
    vector<VectorXd> omega(Nt+1, VectorXd(num_of_FE_nodes_omega));
    vector<VectorXd> u1(Nt+1, VectorXd(num_of_FE_nodes_u));
    vector<VectorXd> u2(Nt+1, VectorXd(num_of_FE_nodes_u));
    vector<VectorXd> p(Nt+1, VectorXd(num_of_FE_nodes_p));

    // initial phi and u (t = 0)
    for (int i = 0; i < num_of_FE_nodes_phi; i++) {
        phi[0][i] = initial_function_phi(M_partition[i][0], M_partition[i][1]);
    }
    string file_name = "../data/phi/phi0.dat";
    result2dat(phi[0], file_name, M_basis_phi, T_basis_phi);

    for (int i = 0; i < num_of_FE_nodes_u; i++) {
        u1[0][i] = initial_function_u1(M_basis_u[i][0], M_basis_u[i][1]);
        u2[0][i] = initial_function_u2(M_basis_u[i][0], M_basis_u[i][1]);
    }
    
    
    // caculate time independent matrix
    cout << "assemble A1." << endl;
    SparseMatrix<double> A1 = assemble_matrix_A(M_partition, T_partition, T_partition, T_partition,
                                                num_of_elements, num_of_FE_nodes_phi, num_of_FE_nodes_phi,
                                                num_of_local_basis_phi, num_of_local_basis_phi,
                                                basis_type_phi, 0, 0, basis_type_phi, 0, 0);
    cout << "assemble A4." << endl;
    SparseMatrix<double> A4 = assemble_matrix_A(
            M_partition, T_partition, T_basis_omega, T_basis_omega,
            num_of_elements, num_of_FE_nodes_omega, num_of_FE_nodes_omega, num_of_local_basis_omega, num_of_local_basis_omega,
            basis_type_omega, 1, 0, basis_type_omega, 1, 0);
    cout << "assemble A5." << endl;
    SparseMatrix<double> A5 = assemble_matrix_A(
            M_partition, T_partition, T_basis_omega, T_basis_omega,
            num_of_elements, num_of_FE_nodes_omega, num_of_FE_nodes_omega, num_of_local_basis_omega, num_of_local_basis_omega,
            basis_type_omega, 0, 1, basis_type_omega, 0, 1);
    cout << "assemble C1." << endl;
    SparseMatrix<double> C1 = assemble_matrix_A(M_partition, T_partition, T_basis_u, T_basis_u,
                                                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u,
                                                num_of_local_basis_u, num_of_local_basis_u,
                                                basis_type_u, 0, 0, basis_type_u, 0, 0);
    cout << "assemble C6." << endl;
    SparseMatrix<double> C6 = assemble_matrix_A(
            M_partition, T_partition, T_basis_u, T_basis_u,
            num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
            basis_type_u, 1, 0, basis_type_u, 1, 0);
    cout << "assemble C7." << endl;
    SparseMatrix<double> C7 = assemble_matrix_A(
            M_partition, T_partition, T_basis_u, T_basis_u,
            num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
            basis_type_u, 0, 1, basis_type_u, 0, 1);
    cout << "assemble C8." << endl;
    SparseMatrix<double> C8 = assemble_matrix_A(M_partition, T_partition, T_basis_p, T_basis_u,
                                                num_of_elements, num_of_FE_nodes_p, num_of_FE_nodes_u,
                                                num_of_local_basis_p, num_of_local_basis_u,
                                                basis_type_p, 0, 0, basis_type_u, 1, 0);
    cout << "assemble D8." << endl;
    SparseMatrix<double> D8 = assemble_matrix_A(
            M_partition, T_partition, T_basis_p, T_basis_u,
            num_of_elements, num_of_FE_nodes_p, num_of_FE_nodes_u, num_of_local_basis_p, num_of_local_basis_u,
            basis_type_p, 0, 0, basis_type_u, 0, 1);
    cout << "assemble G1." << endl;
    SparseMatrix<double> G1 = assemble_matrix_A(M_partition, T_partition, T_basis_p, T_basis_p,
                                                num_of_elements, num_of_FE_nodes_p, num_of_FE_nodes_p,
                                                num_of_local_basis_p, num_of_local_basis_p,
                                                basis_type_p, 0, 0, basis_type_p, 0, 0);
    cout << "assemble G2." << endl;
    SparseMatrix<double> G2 = assemble_matrix_A(
            M_partition, T_partition, T_basis_u, T_basis_p,
            num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_p, num_of_local_basis_u, num_of_local_basis_p,
            basis_type_u, 1, 0, basis_type_p, 0, 0);
    cout << "assemble G3." << endl;
    SparseMatrix<double> G3 = assemble_matrix_A(
            M_partition, T_partition, T_basis_u, T_basis_p,
            num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_p, num_of_local_basis_u, num_of_local_basis_p,
            basis_type_u, 0, 1, basis_type_p, 0, 0);
//    sparse2csv(C1, "C1.csv");
//    sparse2csv(C6, "C6.csv");
//    sparse2csv(C7, "C7.csv");
//    sparse2csv(C8, "C8.csv");


    // iterate every time step
//    Nt = 2;
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - begin);
    std::cout << "before time step: " << elapsed.count() << " milliseconds" << std::endl;
    for (int i = 1; i <= Nt; i++) {
        auto time_begin = std::chrono::steady_clock::now();
        ostringstream stime;
        stime <<  (i*dt);
        string str_time = stime.str();
        cout << "time: "  << str_time << endl;
        // solve phi equation
        cout << "solve phi equation " << endl;
        SparseMatrix<double> A2 = assemble_matrix_A_FE(
                u1[i-1], basis_type_u, 0, 0, M_partition, T_partition, T_basis_phi, T_basis_phi, T_basis_u,
                num_of_elements, num_of_FE_nodes_phi, num_of_FE_nodes_phi, num_of_local_basis_phi, num_of_local_basis_phi,
                basis_type_phi, 1, 0, basis_type_phi, 0, 0);

        SparseMatrix<double> A3 = assemble_matrix_A_FE(
                u2[i-1], basis_type_u, 0, 0, M_partition, T_partition, T_basis_phi, T_basis_phi, T_basis_u,
                num_of_elements, num_of_FE_nodes_phi, num_of_FE_nodes_phi, num_of_local_basis_phi, num_of_local_basis_phi,
                basis_type_phi, 0, 1, basis_type_phi, 0, 0);


        VectorXd b1 = load_vector_FE_b(phi[i-1], basis_type_phi, 0, 0,
                                      M_partition, T_partition, T_basis_phi, T_basis_phi,
                                      num_of_elements, num_of_FE_nodes_phi, num_of_local_basis_phi, basis_type_phi, 0, 0);

        SparseMatrix<double> B1 = assemble_matrix_A_2FE(phi[i-1], basis_type_phi, 0, 0, phi[i-1], basis_type_phi, 0, 0,
                                                        M_partition, T_partition, T_basis_phi, T_basis_phi, T_basis_phi, T_basis_phi,
                                                        num_of_elements, num_of_FE_nodes_phi, num_of_FE_nodes_phi, num_of_local_basis_phi, num_of_local_basis_phi,
                                                        basis_type_phi, 0, 0, basis_type_phi, 0, 0);
        SparseMatrix<double> B2 = A1;
        SparseMatrix<double> B3 = A4;
        SparseMatrix<double> B4 = A5;
        SparseMatrix<double> B5 = A1;

        VectorXd b2(num_of_FE_nodes_omega);
        for (int k = 0; k < num_of_FE_nodes_omega; k++) {
            b2[k] = 0;
        }

        MatrixXd temp_A(num_of_FE_nodes_phi+num_of_FE_nodes_omega, num_of_FE_nodes_phi+num_of_FE_nodes_omega);
        VectorXd right_side_b(num_of_FE_nodes_omega + num_of_FE_nodes_phi);
//        cout << A1 + dt * A2 + dt * A3 << endl;
//        cout << "-------------------------------" << endl;
//        cout << (dt * gamma * A4 + dt * gamma * A5) << endl;
//        cout << "-------------------------------" << endl;
//        cout << (-1/(eta*eta) * (B1 - B2) - B3 - B4) << endl;
//        cout << "-------------------------------" << endl;
//        cout << B5 << endl;
//        cout << "-------------------------------" << endl;

        temp_A << (MatrixXd)(A1 + dt * A2 + dt * A3), (MatrixXd)(dt * gamma * A4 + dt * gamma * A5), (MatrixXd)(-1/(eta*eta) * (B1 - B2) - B3 - B4), (MatrixXd)(B5);
        SparseMatrix<double> left_side_A = temp_A.sparseView();
        left_side_A.makeCompressed();
        cout << "assemble matrix A done." << endl;

//        string file_A_phi = "../data/A_phi/A_phi" + str_time + ".csv";
//        sparse2csv(left_side_A, file_A_phi);

        right_side_b << b1, b2;
//        string file_b_phi = "../data/b_phi/b_phi" + str_time + ".csv";
//        vectro2csv(right_side_b, file_b_phi);

        cout << "load vector b done." << endl;
//        BiCGSTAB<SparseMatrix<double>> solver;
//        LeastSquaresConjugateGradient <SparseMatrix<double>> solver;
        ConjugateGradient <SparseMatrix<double>, Upper> solver;
        solver.compute(left_side_A);
        VectorXd X = solver.solve(right_side_b);

        for (int k = 0; k < num_of_FE_nodes_phi; k++) {
            phi[i][k] = X[k];
//            cout << phi[i][k] << "\t";
//            if ((k+1) % (Nx_partition+1) == 0) {
//                cout << endl;
//            }
        }
        for (int k = 0; k < num_of_FE_nodes_omega; k++) {
            omega[i][k] = X[num_of_FE_nodes_phi+k];
        }

        // solve u & p equation
        cout << "solve u, p equation." << endl;

        SparseMatrix<double> C2 = assemble_matrix_A_FE(
                u1[i-1], basis_type_u, 0, 0, M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
                basis_type_u, 1, 0, basis_type_u, 0, 0);
        cout << "assemble C2." << endl;
        SparseMatrix<double> C3 = assemble_matrix_A_FE(
                u2[i-1], basis_type_u, 0, 0, M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
                basis_type_u, 0, 1, basis_type_u, 0, 0);
        cout << "assemble C3." << endl;
        SparseMatrix<double> C4 = assemble_matrix_A_FE(
                u1[i-1], basis_type_u, 1, 0, M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
                basis_type_u, 0, 0, basis_type_u, 0, 0);
        cout << "assemble C4." << endl;
        SparseMatrix<double> C5 = assemble_matrix_A_FE(
                u2[i-1], basis_type_u, 0, 1, M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
                basis_type_u, 0, 0, basis_type_u, 0, 0);
        cout << "assemble C5." << endl;

        VectorXd e1 = load_vector_2FE_b(omega[i], basis_type_omega, 0, 0,
                                        phi[i], basis_type_phi, 1, 0,
                                       M_partition, T_partition, T_basis_u, T_basis_omega, T_basis_phi,
                                       num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);
        VectorXd e2 = load_vector_2FE_b(omega[i], basis_type_omega, 0, 0,
                                        phi[i], basis_type_phi, 0, 1,
                                        M_partition, T_partition, T_basis_u, T_basis_omega, T_basis_phi,
                                        num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);
        VectorXd e3 = load_vector_2FE_b(u1[i], basis_type_u, 0, 0,
                                        u1[i], basis_type_u, 1, 0,
                                        M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                                        num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);
        VectorXd e4 = load_vector_2FE_b(u2[i], basis_type_u, 0, 0,
                                        u1[i], basis_type_u, 0, 1,
                                        M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                                        num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);
        VectorXd e5 = load_vector_FE_b(u1[i-1], basis_type_u, 0, 0,
                                       M_partition, T_partition, T_basis_u, T_basis_u,
                                       num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);

        SparseMatrix<double> D1 = C1;
        SparseMatrix<double> D2 = C2;
        SparseMatrix<double> D3 = C3;
        SparseMatrix<double> D4 = assemble_matrix_A_FE(
                u2[i-1], basis_type_u, 1, 0, M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
                basis_type_u, 0, 0, basis_type_u, 0, 0);
        cout << "assemble D4." << endl;
        SparseMatrix<double> D5 = assemble_matrix_A_FE(
                u2[i-1], basis_type_u, 0, 1, M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                num_of_elements, num_of_FE_nodes_u, num_of_FE_nodes_u, num_of_local_basis_u, num_of_local_basis_u,
                basis_type_u, 0, 0, basis_type_u, 0, 0);
        cout << "assemble D5." << endl;
        SparseMatrix<double> D6 = C6;
        SparseMatrix<double> D7 = C7;

        VectorXd o1 = e1;
        VectorXd o2 = e2;
        VectorXd o3 = load_vector_2FE_b(u1[i], basis_type_u, 0, 0,
                                        u2[i], basis_type_u, 1, 0,
                                        M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                                        num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);
        VectorXd o4 = load_vector_2FE_b(u2[i], basis_type_u, 0, 0,
                                        u2[i], basis_type_u, 0, 1,
                                        M_partition, T_partition, T_basis_u, T_basis_u, T_basis_u,
                                        num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);
        VectorXd o5 = load_vector_FE_b(u2[i-1], basis_type_u, 0, 0,
                                       M_partition, T_partition, T_basis_u, T_basis_u,
                                       num_of_elements, num_of_FE_nodes_u, num_of_local_basis_u, basis_type_u, 0, 0);


        MatrixXd temp_A2(num_of_FE_nodes_u*2 + num_of_FE_nodes_p, num_of_FE_nodes_u*2 + num_of_FE_nodes_p);
        VectorXd right_side_b2(num_of_FE_nodes_u*2 + num_of_FE_nodes_p);

        SparseMatrix<double> A11 = rho_0 / dt * C1 + rho_0 * C2 + rho_0 * C3 + rho_0 * C5 + mu_0 * C6 + mu_0 * C7;
        SparseMatrix<double> A12 = rho_0 * C4;
        SparseMatrix<double> A13 = C8;
        SparseMatrix<double> A21 = rho_0 * D4;
        SparseMatrix<double> A22 = rho_0/dt * D1 + rho_0 * D2 + rho_0 * D3 + rho_0 * D5 + mu_0 * D6 + mu_0 * D7;
        SparseMatrix<double> A23 = D8;
        SparseMatrix<double> A31 = -1 * G2;
        SparseMatrix<double> A32 = -1 * G3;
        SparseMatrix<double> A33 = G1;
//        sparse2csv(A11, "A11.csv");
//        cout << A11.rows() << " " << A11.cols() << endl;
//        cout << A12.rows() << " " << A12.cols() << endl;
//        cout << A13.rows() << " " << A13.cols() << endl;
//        cout << A21.rows() << " " << A21.cols() << endl;
//        cout << A22.rows() << " " << A22.cols() << endl;
//        cout << A23.rows() << " " << A23.cols() << endl;
//        cout << A31.rows() << " " << A31.cols() << endl;
//        cout << A32.rows() << " " << A32.cols() << endl;
//        cout << A33.rows() << " " << A33.cols() << endl;
//        cout << temp_A2.rows() << " " << temp_A2.cols() << endl;

        temp_A2 << (MatrixXd)(A11), (MatrixXd)(A12), (MatrixXd)(A13),
                (MatrixXd)(A21), (MatrixXd)(A22), (MatrixXd)(A23),
                (MatrixXd)(A31), (MatrixXd)(A32), (MatrixXd)(A33);
        cout << "assemble matrix A done." << endl;
//        temp_A2 << (MatrixXd)(rho_0 / dt * C1 + rho_0 * C2 + rho_0 * C3 + rho_0 * C5 + mu_0 * C6 + mu_0 * C7), (MatrixXd)(rho_0 * C4), (MatrixXd)(C8),
//                   (MatrixXd)(rho_0 * D4), (MatrixXd)(rho_0/dt * D1 + rho_0 * D2 + rho_0 * D3 + rho_0 * D5 + mu_0 * D6 + mu_0 * D7), (MatrixXd)(D8),
//                   (MatrixXd)(-1 * G2), (MatrixXd)(-1 * G3), (MatrixXd)(G1);
        SparseMatrix<double> left_side_A2 = temp_A2.sparseView();
        left_side_A2.makeCompressed();
        VectorXd b21(num_of_FE_nodes_u), b22(num_of_FE_nodes_u), b23(num_of_FE_nodes_p);
        b21 = lambda * e1 + lambda * e2 + rho_0 * e3 + rho_0 * e4 + rho_0 / dt * e5;
        b22 = lambda * o1 + lambda * o2 + rho_0 * o3 + rho_0 * o4 + rho_0 / dt * o5;
        for (int k = 0; k < num_of_FE_nodes_p; k++) {
            b23[k] = 0;
        }
//        cout << "b size:" << endl;
//        cout << b21.size() << endl;
//        cout << b22.size() << endl;
//        cout << b23.size() << endl;
//        cout << right_side_b2.size() << endl;
        right_side_b2 << b21, b22, b23;
        cout << "load vector b done." << endl;

        // write into file
//        string file_A_up = (string)("../data/A_up/A_up") + str_time +  (string)(".csv");
//        cout << file_A_up << endl;
//        sparse2csv(left_side_A2, file_A_up);

        // treat boundary condition
//        cout << "treat boundary: left side" << endl;
//        left_side_A2 = treat_boundary_A_three_matrix(boundary_nodes_u, num_of_FE_nodes_u, boundary_nodes_u, num_of_FE_nodes_u, boundary_nodes_p, left_side_A2);
//        cout << "treat boundary: right side" << endl;
//        right_side_b2 = treat_boundary_b_three(Omega, right_side_b2,
//                               boundary_function_u1, boundary_nodes_u, M_basis_u, num_of_FE_nodes_u,
//                               boundary_function_u2, boundary_nodes_u, M_basis_u, num_of_FE_nodes_u,
//                               boundary_function_p, boundary_nodes_p, M_basis_p);


//        string file_A_boundary = (string)("../data/A_up_boundary/A_up_boundary") + str_time +  (string)(".csv");
//        cout << file_A_boundary<< endl;
//        sparse2csv(left_side_A2, file_A_boundary);

//        string file_b_up_boundary = (string)("../data/b_up_boundary/b_up_boundary") + str_time +  (string)(".csv");
//        cout << file_b_up_boundary << endl;
//        vectro2csv(right_side_b2, file_b_up_boundary);

        // solve the equation
        LeastSquaresConjugateGradient <SparseMatrix<double>> solver2;
//        BiCGSTAB<SparseMatrix<double>> solver2;
        cout << "solver compute left side." << endl;
        solver2.compute(left_side_A2);
        cout << "solver solve with right side" << endl;
        auto solve_begin = std::chrono::steady_clock::now();
        VectorXd X2 = solver2.solve(right_side_b2);
        auto solve_now = std::chrono::steady_clock::now();
        auto solve_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(solve_now - solve_begin);
        std::cout << "time solve: " << solve_elapsed.count() << " milliseconds" << std::endl;

        cout << "get the result of u1,u2,p." << endl;
        for (int k = 0; k < num_of_FE_nodes_u; k++) {
            u1[i][k] = X2[k];
            u2[i][k] = X2[num_of_FE_nodes_u+k];
        }
        for (int k = 0; k < num_of_FE_nodes_p; k++) {
            p[i][k] = X2[num_of_FE_nodes_u * 2 + k];
        }
        if (i <= 10) {
            string file_phi_dat = (string)("../data/phi/phi") + str_time +  (string)(".dat");
            result2dat(phi[i], file_phi_dat, M_basis_phi, T_basis_phi);
//            string file_phi_csv = (string)("../data/phi/phi") + str_time +  (string)(".csv");
//            vectro2csv(phi[i], file_phi_csv);

            string file_omega = (string)("../data/omega/omega_") + str_time +  (string)(".csv");
            vectro2csv(omega[i], file_omega);

            string file_u1 = (string)("../data/u1/u1_") + str_time +  (string)(".dat");
            result2dat(u1[i], file_u1, M_basis_u, T_basis_u);

            string file_u2 = (string)("../data/u2/u2_") + str_time +  (string)(".dat");
            result2dat(u2[i], file_u2, M_basis_u, T_basis_u);
        } else {
            if (i % 10 == 0) {
                string file_phi_dat = (string)("../data/phi/phi") + str_time +  (string)(".dat");
                result2dat(phi[i], file_phi_dat, M_basis_phi, T_basis_phi);

                string file_omega = (string)("../data/omega/omega_") + str_time +  (string)(".csv");
                vectro2csv(omega[i], file_omega);

                string file_u1 = (string)("../data/u1/u1_") + str_time +  (string)(".dat");
                result2dat(u1[i], file_u1, M_basis_u, T_basis_u);

                string file_u2 = (string)("../data/u2/u2_") + str_time +  (string)(".dat");
                result2dat(u2[i], file_u2, M_basis_u, T_basis_u);
            }
        }
        auto time_now = std::chrono::steady_clock::now();
        auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_now - time_begin);
        std::cout << "one time step cost: " << time_elapsed.count() << " milliseconds" << std::endl;
        cout << "---------------------------" << endl;
    }

    return 0;

}