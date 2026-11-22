#include <iostream>
#include <cmath>
#include <tuple>
#include <cassert>
#include <fstream>
#include <iomanip>

#define DEBUG


void bline(){
    std::cout << "\n--------------------------\n\n"; 
}

struct tri_diag{
    std::vector<double> main;
    std::vector<double> upper;
    std::vector<double> lower;

    void print_mat(){
        int width = 8;
    std::cout << std::setw(width) << main[0] << std::setw(width) << upper[0] << std::setw(width) << "." << std::setw(width) << "."  << std::setw(width) << ".\n";
    std::cout << std::setw(width) << lower[0] << std::setw(width) << main[1] << std::setw(width) << upper[1] << std::setw(width) << "." << std::setw(width) << ".\n";
    std::cout << std::setw(width) << "."  <<  std::setw(width) << lower[1] << std::setw(width) << "." << std::setw(width) << "." << std::setw(width) << " .\n";
    std::cout << std::setw(width) << "."  <<  std::setw(width) << "."  << std::setw(width) << "." << std::setw(width) << "." << std::setw(width) << " .\n";
    std::cout << std::setw(width) << "."  <<  std::setw(width) << "."  << std::setw(width) << "." << std::setw(width) << upper[upper.size()-2] << std::setw(width) << " .\n";
    std::cout << std::setw(width) << "." <<  std::setw(width) << "."  << std::setw(width) << lower[lower.size()-2] << std::setw(width) << main[main.size()-2] << std::setw(width) << upper.back() << "\n";
    std::cout << std::setw(width) << "." <<  std::setw(width) << "."  << std::setw(width) << "." << std::setw(width) << lower.back() << std::setw(width) << main.back() << "\n";
    }
};

double piecewiese_linear_1(std::vector<double> grid, double h, int step){
    return (grid[step] - grid[step-1])/h;
}

double piecewiese_linear_2(std::vector<double> grid, double h, int step){
    return (grid[step+1] - grid[step])/h;
}

double ana_righthand_side(double x){
    return -1/2 * std::pow((x - (1/2) ), 2) + 1/8;
}

std::vector<double> create_grid(double h, int N){
    std::vector<double> grid;
    grid.resize(N+1);
    for(int i = 0; i < N+1; i++){
       grid[i] = h*i;
    }
    #ifdef DEBUG
    std::cout << "Number of Grid Points: " << grid.size() << "\n";
    std::cout << "Step Size: " << h << "\n";
    std::cout << "First entry of Grid vec: " << grid.front() << "\n";
    std::cout << "Last entry of Grid vec: " << grid.back() << "\n";
    #endif
    return grid;
}

tri_diag create_mat(double h, int N){
    tri_diag mat;
    mat.main.resize(N-1);
    mat.lower.resize(N-2);
    mat.upper.resize(N-2);
    for(int i = 0; i < N-1; i++){
        mat.main[i] = 2/h;
    }
    for(int i = 0; i < N-2; i++){
        mat.lower[i] = -1/h;
        mat.upper[i] = -1/h;
    }
    #ifdef DEBUG
    bline();
    std::cout << "Size of main diag mat: " << mat.main.size() << " |First Element: " << mat.main[0];
    std::cout << " |Last Element: " << mat.main.back() << "\n";
    std::cout << "Size of lower diag mat: " << mat.lower.size() << " |First Element: " << mat.lower[0] << "\n";
    std::cout << "Size of upper diag mat: " << mat.upper.size() << " |First Element: " << mat.upper[0] << "\n";
    bline();
    #endif
    return mat;
}

std::vector<double> create_b(double h, int N){
    std::vector<double> b(N-1,h);
    #ifdef DEBUG
    std::cout << "Size of b: " << b.size() << " |First Element: " << b[0] << "\n"; 

    #endif
    return b;
}



std::vector<double> thomas_alg(tri_diag mat, std::vector<double> b){
    std::vector<double> U(mat.main.size(), 0);
    double m;

    for(int i = 1; i < mat.main.size(); i++){
        m = mat.lower[i-1] / mat.main[i-1];
        mat.main[i] = mat.main[i] - m * mat.upper[i-1];
        b[i] = b[i] - m * b[i-1];
        mat.lower[i-1] = 0;
    }
    
    U.back() = b.back() / mat.main.back();

    #ifdef DEBUG
    bline();
    mat.print_mat();
    bline();
    std::cout << "Size of result vec: " << U.size() << " |Last entry after forwad calc: ";
    std::cout << U.back() << "\n";
    std::cout << "First element of lower: " << mat.lower[0];
    std::cout << " |Last element of lower: " << mat.lower.back() << "\n";
    #endif

    U.back() = b.back()/mat.main.back();
    for(int i = mat.main.size()-1; i > 0; i--){
        U[i-1] = (b[i-1] - mat.upper[i-1] * U[i]) / mat.main[i-1]; 
    }
    return U;
}


int main(int argc, char* argv[]){
    assert(argc == 2);
    int step_size = std::stoi(argv[1]);
    assert(step_size < 6);
    std::vector<double> h_vec {0.1,0.01,0.001,0.0001,0.00001,0.000001};
    int N = 1/h_vec[step_size];

    tri_diag matrices = create_mat(h_vec[step_size], N);   
    #ifdef DEBUG
    std::cout << "Tri Diag mat: \n";
    matrices.print_mat();
    bline();
    std::cout << "Step sizes in step size vec: " << "\n";
    for(size_t i = 0; i < h_vec.size(); i++){
        std::cout << i+1 << ". " << h_vec[i] <<"\n";
    }
    bline();
    #endif
    std::vector<double> grid1 = create_grid(h_vec[step_size], N);

    
    std::vector<double> b1 = create_b(h_vec[step_size], N);
    std::vector<double> result1 = thomas_alg(matrices, b1);
    #ifdef DEBUG
    bline();
    std::cout << "First element of result: " << result1[0] << " |Centered element of result: " << result1[static_cast<int>(result1.size()/2)];
    std::cout << " |Last element of result: " << result1.back() << "\n";
    #endif
    std::ofstream file1("solution1.dat");

    for (int i = 0; i < N - 1; ++i) {
        double x = (i + 1) * h_vec[step_size];
        file1 << x << " " << result1[i] << "\n";
    }

    file1.close();

    std::vector<double> b2(N,0);
    double temp{};
    int middle{};
    for(int i = 0; i < N; i++){
        if(std::abs(temp-0.5) <1e-10){
            b2[i] = 2;
            middle = i;
        }
        temp = temp + h_vec[step_size];
    }

    #ifdef DEBUG
    bline();
    std::cout << "Tri Diag for b: " << "\n";
    matrices.print_mat();
    
    bline();
    std::cout << "Shape of b for Dirac dist: " << "\n";
    std::cout << b2[0] << "\n";
    std::cout << b2[1] << "\n";
    std::cout << "." << "\n";
    std::cout << b2[middle-1] << "\n";
    std::cout << b2[middle] << "\n";
    std::cout << b2[middle+1] << "\n";
    std::cout << "." << "\n";
    std::cout << b2[b2.size()-2] << "\n";
    std::cout << b2.back() << "\n";

    #endif
    std::vector<double> result2 = thomas_alg(matrices, b2);

    #ifdef DEBUG
    bline();
    std::cout << "First element of result: " << result2[0] << " |Centered element of result: " << result2[static_cast<int>(result1.size()/2)];
    std::cout << " |Last element of result: " << result2.back() << "\n";
    #endif
    std::ofstream file2("solution2.dat");

    for (int i = 0; i < N - 1; ++i) {
        double x = (i + 1) * h_vec[step_size];
        file2 << x << " " << result2[i] << "\n";
    }

    file2.close();



    return 0;
}
