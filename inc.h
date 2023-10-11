#include <string>

double f(int s, int n, int i, int j);
void input_matrix(int s, int n, std::vector<double>* matrix);
void input_b(int n, const std::vector<double>& matrix, std::vector<double>* b);
bool read_file(const std::string& name_file, int n, std::vector<double>* matrix);
void output(int n, int r, int l, const std::vector<double>& vec);
bool solution(int n, int m, std::vector<double>* matrix, std::vector<double>* b,
    std::vector<double>* x, std::vector<int>* block_rows, std::vector<double>* block1, std::vector<double>* block2);
double matrix_norm(int n, int m, const std::vector<double>& matrix);
void matr_prod(int n, int m, int k, const std::vector<double>& a, 
    const std::vector<double>& b, std::vector<double>* c);
void subtract_vector(int n, const std::vector<double>& a, const std::vector<double>& b, std::vector<double>* c);
double vector_norm(const std::vector<double>& vec);
void get_block(int i, int j, int n, int m, int k, int l,
 const std::vector<double>& matrix, std::vector<double>* block1);
void put_block(int i, int j, int n, int m, int k, int l, const std::vector<double>& block, std::vector<double>* matrix);
double r1_eval(const std::vector<double>& matrix, const std::vector<double>& x, const std::vector<double>& b, std::vector<double>* c, std::vector<double>* d);
double r2_eval(const std::vector<double>& x);
bool is_inv(int m, std::vector<double>* matrix, double a_norm);
std::vector<int> inverse_matrix(int m, std::vector<double>* matrix, std::vector<double>* identity);
void subtract_matrix_inplace(int n, int m, std::vector<double>* a, const std::vector<double>& b);
std::vector<double> get_vector(int i, int m, int k, int l, const std::vector<double>& b);
void put_vector(int i, int m, int k, int l, const std::vector<double>& b_i, std::vector<double>* b);
std::vector<double> matrix_product(int n, int m, int k, const std::vector<double>& a,
    const std::vector<double>& b, const std::vector<int>& inv_rows);
    