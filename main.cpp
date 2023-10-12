#include <iostream>
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

    double* matrix = new double[n*n];
    double* b = new double[n];

    if (s) {
        input_matrix(s, n, matrix);
    } else {
        std::string name_file = argv[argc - 1];
        if (!read_file(name_file, n, matrix)) {
            std::cout << "Проблемы с чтением файла" << "\n";

            delete[] matrix;
            delete[] b;
            return -1;
        }
    }

    input_b(n, matrix, b);

    double* x = new double[n];
    double* block1 = new double[m*m];
    double* block2 = new double[m*m];
    double* block3 = new double[m*m];
    int* rows = new int[m];

    int k = n / m;
    int l = n % m;
    int h = l ? k + 1 : k; 
    int* block_rows = new int[h];

    auto start = high_resolution_clock::now();
    bool ok = solution(n, m, matrix, b, x, block_rows, rows, block1, block2, block3);
    auto stop = high_resolution_clock::now();
    duration<double> diff = stop - start;
    double t1 = diff.count();

    if (!ok) {
        printf (
            "%s : Task = %d Res1 = %d Res2 = %d T1 = %.2f T2 = %d S = %d N = %d M = %d\n",
            argv[0], task, -1, -1, t1, 0, s, n, m
        );    

        delete[] matrix;
        delete[] b;
        delete[] x;
        delete[] block1;
        delete[] block2;
        delete[] block3;
        delete[] rows;
        delete[] block_rows;
        return -1;
    } 

    if (s) {
        input_matrix(s, n, matrix);
    } else {
        std::string name_file = argv[argc - 1];
        if (!read_file(name_file, n, matrix)) {
            std::cout << "Проблемы с чтением файла" << "\n";

            delete[] matrix;
            delete[] b;
            delete[] x;
            delete[] block1;
            delete[] block2;
            delete[] block3;
            delete[] rows;
            delete[] block_rows;
            return -1;
        }
    }

    input_b(n, matrix, b);   

    double* c = new double[n];
    double* d = new double[n];
    start = high_resolution_clock::now();
    double r1 = r1_eval(n, matrix, x, b, c, d);
    double r2 = r2_eval(n, x);
    stop = high_resolution_clock::now();
    diff = stop - start;
    double t2 = diff.count();

    std::cout << "A:\n";
    output(n, r, n, matrix);
    std::cout << "\n";

    std::cout << "x:\n";
    output(n, r, 1, x);
    std::cout << "\n";

    printf (
        "%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n",
        argv[0], task, r1, r2, t1, t2, s, n, m
    );


    delete[] matrix;
    delete[] b;
    delete[] x;
    delete[] block1;
    delete[] block2;
    delete[] block3;
    delete[] rows;
    delete[] block_rows;
    delete[] c;
    delete[] d;
    return 0;
}
