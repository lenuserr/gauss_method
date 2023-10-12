#include <string>

double f(int s, int n, int i, int j);
void input_matrix(int s, int n, double* matrix);
void input_b(int n, double* matrix, double* b);
bool read_file(const std::string& name_file, int n, double* matrix);
void output(int n, int r, int l, double* vec);
double r1_eval(int n, double* matrix, double* x, double* b, double* c, double* d);
double r2_eval(int n, double* x);
double matrix_norm(int n, int m, double* matrix);
void matrix_product(int n, int m, int k, double* a, double* b, double* res, int* inv_rows);
void matr_prod(int n, int m, int k, double* a, double* b, double* res);
void subtract_vector(int n, double* a, double* b, double* c);
double vector_norm(int n, double* vec);
void get_block(int i, int j, int n, int m, int k, int l, double* matrix, double* block1);
void put_block(int i, int j, int n, int m, int k, int l, double* block, double* matrix);
void get_vector(int i, int m, int k, int l, double* b, double* block3);
void put_vector(int i, int m, int k, int l, double* b_i, double* b);
void subtract_matrix_inplace(int n, int m, double* a, double* b);
bool inverse_matrix(int m, double* matrix, double* identity, int* rows, double a_norm);
bool is_inv(int m, double* matrix, double a_norm);
bool solution(int n, int m, double* matrix, double* b, double* x, 
    int* block_rows, int* rows, double* block1, double* block2, double* block3);