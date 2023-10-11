#include <iostream>
#include <vector>
#include <chrono>
#include <cstdio>
#include "inc.h"


using namespace std::chrono;

int main(int argc, char* argv[]) {
    const int task = 9;

    int n = std::stoi(argv[1]);
    int m = std::stoi(argv[2]);
    int r = std::stoi(argv[3]);
    int s = std::stoi(argv[4]);

    std::vector<double> matrix(n * n);
    std::vector<double> b(n);

    if (s) {
        input_matrix(s, n, &matrix);
    } else {
        std::string name_file = argv[argc - 1];
        if (!read_file(name_file, n, &matrix)) {
            std::cout << "Проблемы с чтением файла" << "\n";
            return -1;
        }
    }

    input_b(n, matrix, &b);

    auto copy_matrix = matrix;
    auto copy_b = b;

    std::vector<double> x(n);

    auto start = high_resolution_clock::now();
    bool ok = solution(n, m, &matrix, &b, &x);
    auto stop = high_resolution_clock::now();
    duration<double> diff = stop - start;
    double t1 = diff.count();
    
    std::cout << "A:\n";
    output(n, r, n, copy_matrix);
    std::cout << "\n";

    if (!ok) {
        printf (
            "%s : Task = %d Res1 = %d Res2 = %d T1 = %.2f T2 = %d S = %d N = %d M = %d\n",
            argv[0], task, -1, -1, t1, 0, s, n, m
        );    

        return -1;
    }

    std::cout << "x:\n";
    output(n, r, 1, x);
    std::cout << "\n";    

    std::vector<double> c(n);
    std::vector<double> d(n);
    start = high_resolution_clock::now();
    double r1 = r1_eval(copy_matrix, x, copy_b, &c, &d);
    double r2 = r2_eval(x);
    stop = high_resolution_clock::now();
    diff = stop - start;
    double t2 = diff.count();

    printf (
        "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
        argv[0], task, r1, r2, t1, t2, s, n, m
    );
    
    return 0;
}
