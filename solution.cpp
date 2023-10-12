#include <cmath>
#include "inc.h"
#define EPS 1e-14

bool solution(int n, int m, double* matrix, double* b, double* x, 
    int* block_rows, int* rows, double* block1, double* block2, double* block3) {
    
    int k = n / m;
    int l = n % m;
    int h = l ? k + 1 : k; 
    for (int i = 0; i < h; ++i) {
        block_rows[i] = i;
    }

    double a_norm = matrix_norm(n, n, matrix);
    for (int i = 0; i < k; ++i) {
        double max_norm_block = -1;
        int row_max_block = -1;
        for (int j = i; j < k; ++j) {
            get_block(block_rows[j], i, n, m, k, l, matrix, block1);
            double block_norm = matrix_norm(m, m, block1);
            if (block_norm - EPS * a_norm > max_norm_block && is_inv(m, block1, a_norm, rows)) {
                max_norm_block = block_norm;
                row_max_block = j;
            }
        }

        if (row_max_block == -1) {
            return false;
        }

        std::swap(block_rows[i], block_rows[row_max_block]);
        get_block(block_rows[i], i, n, m, k, l, matrix, block1); 

        inverse_matrix(m, block1, block2, rows, a_norm); 

        for (int s = i + 1; s < k; ++s) { 
            get_block(block_rows[i], s, n, m, k, l, matrix, block1); 
            matrix_product(m, m, m, block2, block1, block3, rows); 
            put_block(block_rows[i], s, n, m, k, l, block3, matrix);
        }

        if (l) {
            get_block(block_rows[i], k, n, m, k, l, matrix, block1); 
            matrix_product(m, m, l, block2, block1, block3, rows);
            put_block(block_rows[i], k, n, m, k, l, block3, matrix);
        }

        matrix_product(m, m, 1, block2, b + m * block_rows[i], block3, rows); 
        put_vector(block_rows[i], m, k, l, block3, b);

        for (int q = i + 1; q < h; ++q) {
            get_block(block_rows[q], i, n, m, k, l, matrix, block1); 
            int multiplier_rows = q < k ? m : l;
            for (int r = i + 1; r < h; ++r) { 
                get_block(block_rows[i], r, n, m, k, l, matrix, block2); 
                int block_cols = r < k ? m : l;
                matr_prod(multiplier_rows, m, block_cols, block1, block2, block3);
                get_block(block_rows[q], r, n, m, k, l, matrix, block2);
                subtract_matrix_inplace(multiplier_rows, block_cols, block2, block3);
                put_block(block_rows[q], r, n, m, k, l, block2, matrix);
            }

            matr_prod(multiplier_rows, m, 1, block1, b + m*block_rows[i], block3);
            subtract_matrix_inplace(1, multiplier_rows, b + m*block_rows[q], block3);
        }
    }

    if (l) {
        get_block(block_rows[k], k, n, m, k, l, matrix, block1);  

        if (!inverse_matrix(l, block1, block2, rows, a_norm)) {
            return false;
        } 

        matrix_product(l, l, 1, block2, b + m*block_rows[k], block3, rows);
        put_vector(k, m, k, l, block3, x);
    }

    for (int i = k - 1; i >= 0; --i) {
        for (int vv = 0; vv < m; ++vv) { 
            block2[vv] = 0;
        }

        for (int j = i + 1; j < h; ++j) {
            get_block(block_rows[i], j, n, m, k, l, matrix, block1); 
            int cols = j < k ? m : l;
            matr_prod(m, cols, 1, block1, x + m*j, block3);

            for (int p = 0; p < m; ++p) {
                block2[p] += block3[p];
            }
        }

        for (int p = 0; p < m; ++p) {
            b[(m*block_rows[i]) + p] -= block2[p]; 
        }

        put_vector(i, m, k, l, b + m*block_rows[i], x);
    }     

    return true;
}

double matrix_norm(int n, int m, double* matrix) {
    double norm = -1;
    for (int j = 0; j < m; ++j) {
        double sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += std::fabs(matrix[m*i + j]);           
        }

        norm = std::max(norm, sum);
    }

    return norm;
}

void matrix_product(int n, int m, int k, double* a, double* b, double* c, int* inv_rows) {
    
    for (int i = 0; i < n*k; ++i) {
        c[i] = 0;
    }
    
    double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;
    
    int v3 = n%3;
    int h3 = k%3;
    
    for (int i = 0; i < v3; ++i) {
        for (int j = 0; j < h3; ++j) {
            sum00 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*inv_rows[i] + p] * b[k*p + j];
            }
            
            c[k*i + j] = sum00;
        }
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*inv_rows[i] + p] * b[k*p + j];
                sum01 += a[m*inv_rows[i] + p] * b[k*p + j + 1];
                sum02 += a[m*inv_rows[i] + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01;
            c[k*i + j + 2] = sum02;
        }
    }
    
    for (int i = v3; i < n; i+=3) {   
        for (int j = 0; j < h3; ++j) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*inv_rows[i] + p] * b[k*p + j];
                sum01 += a[m*inv_rows[i + 1] + p] * b[k*p + j];
                sum02 += a[m*inv_rows[i + 2] + p] * b[k*p + j];
            }
            c[k*i + j] = sum00;
            c[k*(i + 1) + j] = sum01;
            c[k*(i + 2) + j] = sum02;
        }
        
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0; sum10 = 0; sum11 = 0; sum12 = 0;
            sum20 = 0; sum21 = 0; sum22 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*inv_rows[i] + p] * b[k*p + j];
                sum01 += a[m*inv_rows[i] + p] * b[k*p + j + 1];
                sum02 += a[m*inv_rows[i] + p] * b[k*p + j + 2];
                sum10 += a[m*inv_rows[i + 1] + p] * b[k*p + j];
                sum11 += a[m*inv_rows[i + 1] + p] * b[k*p + j + 1];
                sum12 += a[m*inv_rows[i + 1] + p] * b[k*p + j + 2];
                sum20 += a[m*inv_rows[i + 2] + p] * b[k*p + j];
                sum21 += a[m*inv_rows[i + 2] + p] * b[k*p + j + 1];
                sum22 += a[m*inv_rows[i + 2] + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01; 
            c[k*i + j + 2] = sum02;
            c[k*(i + 1) + j] = sum10;
            c[k*(i + 1)  + j + 1] = sum11;
            c[k*(i + 1) + j + 2] = sum12;
            c[k*(i + 2) + j] = sum20;
            c[k*(i + 2) + j + 1] = sum21;
            c[k*(i + 2) + j + 2] = sum22;
        }
    }    
}

void matr_prod(int n, int m, int k, double* a, double* b, double* c) {
    
    for (int i = 0; i < n*k; ++i) {
        c[i] = 0;
    }
    
    double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;
    
    int v3 = n%3;
    int h3 = k%3;
    
    for (int i = 0; i < v3; ++i) {
        for (int j = 0; j < h3; ++j) {
            sum00 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
            }
            
            c[k*i + j] = sum00;
        }
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*i + p] * b[k*p + j + 1];
                sum02 += a[m*i + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01;
            c[k*i + j + 2] = sum02;
        }
    }
    
    for (int i = v3; i < n; i+=3) {   
        for (int j = 0; j < h3; ++j) {
            sum00 = 0; sum01 = 0; sum02 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*(i + 1) + p] * b[k*p + j];
                sum02 += a[m*(i + 2) + p] * b[k*p + j];
            }
            c[k*i + j] = sum00;
            c[k*(i + 1) + j] = sum01;
            c[k*(i + 2) + j] = sum02;
        }
        
        for (int j = h3; j < k; j+=3) {
            sum00 = 0; sum01 = 0; sum02 = 0; sum10 = 0; sum11 = 0; sum12 = 0;
            sum20 = 0; sum21 = 0; sum22 = 0;
            for (int p = 0; p < m; ++p) {
                sum00 += a[m*i + p] * b[k*p + j];
                sum01 += a[m*i + p] * b[k*p + j + 1];
                sum02 += a[m*i + p] * b[k*p + j + 2];
                sum10 += a[m*(i + 1) + p] * b[k*p + j];
                sum11 += a[m*(i + 1) + p] * b[k*p + j + 1];
                sum12 += a[m*(i + 1) + p] * b[k*p + j + 2];
                sum20 += a[m*(i + 2) + p] * b[k*p + j];
                sum21 += a[m*(i + 2) + p] * b[k*p + j + 1];
                sum22 += a[m*(i + 2) + p] * b[k*p + j + 2];
            }
            c[k*i + j] = sum00;
            c[k*i + j + 1] = sum01; 
            c[k*i + j + 2] = sum02;
            c[k*(i + 1) + j] = sum10;
            c[k*(i + 1)  + j + 1] = sum11;
            c[k*(i + 1) + j + 2] = sum12;
            c[k*(i + 2) + j] = sum20;
            c[k*(i + 2) + j + 1] = sum21;
            c[k*(i + 2) + j + 2] = sum22;
        }
    }    
}

void subtract_vector(int n, double* a, double* b, double* c) {   
    for (int i = 0; i < n; ++i) {
        c[i] = a[i] - b[i];
    }
}

double vector_norm(int n, double* vec) {
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        norm += std::fabs(vec[i]);        
    }

    return norm;
}

void get_block(int i, int j, int n, int m, int k, int l, double* matrix, double* block1) {
    int h = i < k ? m : l;
    int w = j < k ? m : l;
    
    int ind = 0;

    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            block1[ind] = matrix[n * (m * i + p) + m * j + q];
            ind++;
        }
    }
}

void put_block(int i, int j, int n, int m, int k, int l, double* block, double* matrix) {
    int h = i < k ? m : l;
    int w = j < k ? m : l;

    int ind = 0;
    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            matrix[n * (m * i + p) + m * j + q] = block[ind];
            ind++;
        }
    }
}

void put_vector(int i, int m, int k, int l, double* b_i, double* b) {
    int length = i < k ? m : l;
    for (int p = 0; p < length; ++p) {
        b[m*i + p] = b_i[p];
    }
}

void subtract_matrix_inplace(int n, int m, double* a, double* b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            a[m*i + j] -= b[m*i + j];           
        }
    }
}

bool inverse_matrix(int m, double* matrix, double* identity, int* rows, double a_norm) {

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            identity[i * m + j] = (i != j) ? 0 : 1;
        }
    }

    for (int k = 0; k < m; ++k) {
        rows[k] = k;
    }

    for (int i = 0; i < m; ++i) {
        double max_elem = matrix[rows[i] * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs(matrix[rows[j] * m + i]) > std::fabs(max_elem)) {
                max_elem = matrix[rows[j] * m + i];
                row_max_elem = j;
            }
        }

        std::swap(rows[i], rows[row_max_elem]);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        double factor = 1 / max_elem;
        for (int s = 0; s < i; ++s) {
            identity[rows[i] * m + s] *= factor;
        }

        for (int s = i; s < m; ++s) {
            matrix[rows[i] * m + s] *= factor;
            identity[rows[i] * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -matrix[rows[k] * m + i];
            for (int p = 0; p < i + 1; ++p) {
                identity[rows[k] * m + p] += identity[rows[i] * m + p] * multiplier;
            }

            for (int p = i + 1; p < m; ++p) { 
                matrix[rows[k] * m + p] += matrix[rows[i] * m + p] * multiplier;
                identity[rows[k] * m + p] += identity[rows[i] * m + p] * multiplier;
            }
        }
    }

    for (int i = m - 1; i > 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            double multiplier = -matrix[rows[k] * m + i];
            for (int p = 0; p < m; ++p) { 
                identity[rows[k] * m + p] += identity[rows[i] * m + p] * multiplier;
            }
        }
    }

    return true;
}

bool is_inv(int m, double* matrix, double a_norm, int* rows) {
    for (int k = 0; k < m; ++k) {
        rows[k] = k;
    }

    for (int i = 0; i < m; ++i) {
        double max_elem = matrix[rows[i] * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs(matrix[rows[j] * m + i]) - EPS * a_norm > std::fabs(max_elem)) {
                max_elem = matrix[rows[j] * m + i];
                row_max_elem = j;
            }
        }

        std::swap(rows[i], rows[row_max_elem]);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        double factor = 1 / max_elem;
        for (int s = i; s < m; ++s) {
            matrix[rows[i] * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -matrix[rows[k] * m + i];
            for (int p = i + 1; p < m; ++p) { 
                matrix[rows[k] * m + p] += matrix[rows[i] * m + p] * multiplier;
            }
        }
    }

    return true;
}
