#include <fstream>
#include <cmath>

double f(int s, int n, int i, int j) { 
    switch(s) {
        case 1:
            return n - std::max(i, j);
        case 2:
            return std::max(i, j) + 1;
        case 3:
            return std::fabs(i-j);
        case 4:
            return 1. / (i + j + 1);
        default:
            return 1;
    }
}

void input_matrix(int s, int n, double* matrix) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[n * i + j] = f(s, n, i, j);
        }
    }       
}

void input_b(int n, double* matrix, double* b) {
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int k = 0; k <= (n-1) / 2; ++k) {
            sum += matrix[n * i + 2*k];
        }
        b[i] = sum; 
    } 
}

bool read_file(const std::string& name_file, int n, double* matrix) {
    std::ifstream fin(name_file);
    
    if (!fin.is_open()) {
        return false;
    }
    
    double x;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(fin >> x)) {
                return false;  
            }
            
            matrix[n * i + j] = x;
        }
    }
    return true;
}
