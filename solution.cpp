#include <vector>
#include <cmath>
#include "inc.h"
#define EPS 1e-14

// 18:00. Переписываю правильно код и уже в конце меняю просто std::vector<double> на double* и всё.
// 22:00. Продолжаю заниматься этим же. Был перерыв на баги + на еду.

bool solution(int n, int m, std::vector<double>* matrix, std::vector<double>* b, std::vector<double>* x,
 std::vector<int>* block_rows, std::vector<int>* rows, 
 std::vector<double>* block1, std::vector<double>* block2, std::vector<double>* block3) {
    
    int k = n / m;
    int l = n % m;
    int h = l ? k + 1 : k; // h блочных строк у меня.
    for (int i = 0; i < h; ++i) {
        (*block_rows)[i] = i;
    }

    double a_norm = matrix_norm(n, n, *matrix);
    for (int i = 0; i < k; ++i) {
        double max_norm_block = -1;
        int row_max_block = -1;
        for (int j = i; j < k; ++j) {
            get_block((*block_rows)[j], i, n, m, k, l, *matrix, block1);
            double block_norm = matrix_norm(m, m, *block1);
            if (block_norm - EPS * a_norm > max_norm_block && is_inv(m, block1, a_norm)) {
                max_norm_block = block_norm;
                row_max_block = j;
            }
        }

        if (row_max_block == -1) {
            return false;
        }

        std::swap((*block_rows)[i], (*block_rows)[row_max_block]);
        get_block((*block_rows)[i], i, n, m, k, l, *matrix, block1); 

        inverse_matrix(m, block1, block2, rows); 

        // block2 занят под обратную матрицу.

        // умножаем блочную строку слева на обратную матрицу.
        for (int s = i + 1; s < k; ++s) { 
            get_block((*block_rows)[i], s, n, m, k, l, *matrix, block1); 
            matrix_product(m, m, m, *block2, *block1, block3, *rows); 
            put_block((*block_rows)[i], s, n, m, k, l, *block3, matrix);
        }

        if (l) {
            get_block((*block_rows)[i], k, n, m, k, l, *matrix, block1); 
            matrix_product(m, m, l, *block2, *block1, block3, *rows);
            put_block((*block_rows)[i], k, n, m, k, l, *block3, matrix);
        }

        // также не забываем слева умножить нужную часть вектора b слева на inv_max_block.
        get_vector((*block_rows)[i], m, k, l, *b, block1); 
        matrix_product(m, m, 1, *block2, *block1, block3, *rows); 
        put_vector((*block_rows)[i], m, k, l, *block3, b);

        for (int q = i + 1; q < h; ++q) {
            get_block((*block_rows)[q], i, n, m, k, l, *matrix, block1); // multiplier = block1
            int multiplier_rows = q < k ? m : l;
            // multiplier shape: [multiplier_rows, m]
            for (int r = i + 1; r < h; ++r) { 
                get_block((*block_rows)[i], r, n, m, k, l, *matrix, block2); 
                int block_cols = r < k ? m : l;
                // multiplier shape: [multiplier_rows, m]; block shape: [m, block_cols]

                std::vector<double> result(multiplier_rows * block_cols);
                matr_prod(multiplier_rows, m, block_cols, *block1, *block2, &result);
                get_block((*block_rows)[q], r, n, m, k, l, *matrix, block2);
                subtract_matrix_inplace(multiplier_rows, block_cols, block2, result);
                put_block((*block_rows)[q], r, n, m, k, l, *block2, matrix);
            }

            // Аналогичные формулы нужно выполнить для вектора b:
            get_vector((*block_rows)[i], m, k, l, *b, block2);
            std::vector<double> v_i(multiplier_rows);
            matr_prod(multiplier_rows, m, 1, *block1, *block2, &v_i);
            get_vector((*block_rows)[q], m, k, l, *b, block1);
            subtract_matrix_inplace(1, multiplier_rows, block1, v_i);
            put_vector((*block_rows)[q], m, k, l, *block1, b);
        }
    }

    // Осталось написать обратный ход метода Гаусса и всё.
    // костыль для самого маленького блока l x l. я про него забыл, поэтому не хочу весь код выше менять.
    if (l) {
        get_block((*block_rows)[k], k, n, m, k, l, *matrix, block1);  

        auto copy_last_block = *block1;
        if (!is_inv(l, &copy_last_block, a_norm)) {
            return false;
        }
        // потом вот этот copy last block уберу нахой, этот код не нужен.
        // МОЖНО УБРАТЬ, ЕСЛИ INVERSE_MATRIX СДЕЛАТЬ BOOL И ТАКУЮ ЛОГИКУ КАК В IS_INV НА RETURN FALSE.

        inverse_matrix(l, block1, block2, rows); 
        get_vector((*block_rows)[k], m, k, l, *b, block1);
        matrix_product(l, l, 1, *block2, *block1, block3, *rows);
        put_vector(k, m, k, l, *block3, x);
    }

    // Обратный ход метода Гаусса.
    for (int i = k - 1; i >= 0; --i) {
        std::vector<double> res(m);
        for (int j = i + 1; j < h; ++j) {
            get_block((*block_rows)[i], j, n, m, k, l, *matrix, block1); 
            output(m, m, m, *block1);
            int cols = j < k ? m : l;
            std::vector<double> prod(m);
            get_vector(j, m, k, l, *x, block2);  
            matr_prod(m, cols, 1, *block1, *block2, &prod);

            for (int p = 0; p < m; ++p) {
                res[p] += prod[p];
            }
        }

        get_vector((*block_rows)[i], m, k, l, *b, block1); 
        for (int p = 0; p < m; ++p) {
            (*block1)[p] -= res[p]; 
        }

        put_vector(i, m, k, l, *block1, x);
    }     

    return true;
}

double matrix_norm(int n, int m, const std::vector<double>& matrix) {
    // matrix shape: n x m
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

void matrix_product(int n, int m, int k, const std::vector<double>& a, const std::vector<double>& b, 
std::vector<double>* res, const std::vector<int>& inv_rows) {
    int i, j, l;
    double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;

    int v = n % 3;
    int h = k % 3;
    for(i = 0; i < n * k; i++)
          (*res)[i] = 0;
    for(i = 0; i < v; i++) {
            for(j = 0; j < h; j++) {
                    sum00 = 0;
                    for(l = 0; l < m; l++)
                            sum00 += a[inv_rows[i]*m+l] * b[l*k+j];
                    (*res)[i*k+j] += sum00;
            }
            for(; j < k; j += 3) {
                    sum00 = sum01 = sum02 = 0;
                            for(l = 0; l < m; l++) {
                                    sum00 += a[inv_rows[i]*m+l] * b[l*k+j];
                                    sum01 += a[inv_rows[i]*m+l] * b[l*k+j+1];
                                    sum02 += a[inv_rows[i]*m+l] * b[l*k+j+2];
                            }
                    (*res)[i*k+j] += sum00;
                    (*res)[i*k+j+1] += sum01;
                    (*res)[i*k+j+2] += sum02;
            }
    }
    for(; i < n; i += 3) {
            for(j = 0; j < h; j++) {
                    sum00 = sum10 = sum20 = 0;
                    for(l = 0; l < m; l++) {
                            sum00 += a[inv_rows[i]*m+l] * b[l*k+j];
                            sum10 += a[inv_rows[i + 1]*m+l] * b[l*k+j];
                            sum20 += a[inv_rows[i + 2]*m+l] * b[l*k+j];
                    }
                    (*res)[i*k+j] += sum00;
                    (*res)[(i+1)*k+j] += sum10;
                    (*res)[(i+2)*k+j] += sum20;
            }

            for(; j < k; j += 3) {
                    sum00 = sum01 = sum02 = sum10 = sum11 = sum12 = sum20 = sum21 = sum22 = 0;
	for(l = 0; l < m; l++) {
                            sum00 += a[inv_rows[i]*m+l] * b[l*k+j];
                            sum01 += a[inv_rows[i]*m+l] * b[l*k+j+1];
                            sum02 += a[inv_rows[i]*m+l] * b[l*k+j+2];
                            sum10 += a[inv_rows[i + 1]*m+l] * b[l*k+j];
                            sum11 += a[inv_rows[i + 1]*m+l] * b[l*k+j+1];
                            sum12 += a[inv_rows[i + 1]*m+l] * b[l*k+j+2];
                            sum20 += a[inv_rows[i + 2]*m+l] * b[l*k+j];
                            sum21 += a[inv_rows[i + 2]*m+l] * b[l*k+j+1];
                            sum22 += a[inv_rows[i + 2]*m+l] * b[l*k+j+2];
                    }
                    (*res)[i*k+j] += sum00;
                    (*res)[i*k+j+1] += sum01;
                    (*res)[i*k+j+2] += sum02;
                    (*res)[(i+1)*k+j] += sum10;
                    (*res)[(i+1)*k+j+1] += sum11;
                    (*res)[(i+1)*k+j+2] += sum12;
                    (*res)[(i+2)*k+j] += sum20;
                    (*res)[(i+2)*k+j+1] += sum21;
                    (*res)[(i+2)*k+j+2] += sum22; 
            }
    }
}

void matr_prod(int n, int m, int k, const std::vector<double>& a, const std::vector<double>& b, std::vector<double>* res) {
        int i, j, l;
        double sum00, sum01, sum02, sum10, sum11, sum12, sum20, sum21, sum22;

        int v = n % 3;
        int h = k % 3;
        for(i = 0; i < n * k; i++)
              (*res)[i] = 0;
        for(i = 0; i < v; i++) {
                for(j = 0; j < h; j++) {
                        sum00 = 0;
                        for(l = 0; l < m; l++)
                                sum00 += a[i*m+l] * b[l*k+j];
                        (*res)[i*k+j] += sum00;
                }
                for(; j < k; j += 3) {
                        sum00 = sum01 = sum02 = 0;
                                for(l = 0; l < m; l++) {
                                        sum00 += a[i*m+l] * b[l*k+j];
                                        sum01 += a[i*m+l] * b[l*k+j+1];
                                        sum02 += a[i*m+l] * b[l*k+j+2];
                                }
                        (*res)[i*k+j] += sum00;
                        (*res)[i*k+j+1] += sum01;
                        (*res)[i*k+j+2] += sum02;
                }
        }
        for(; i < n; i += 3) {
                for(j = 0; j < h; j++) {
                        sum00 = sum10 = sum20 = 0;
                        for(l = 0; l < m; l++) {
                                sum00 += a[i*m+l] * b[l*k+j];
                                sum10 += a[(i+1)*m+l] * b[l*k+j];
                                sum20 += a[(i+2)*m+l] * b[l*k+j];
                        }
                        (*res)[i*k+j] += sum00;
                        (*res)[(i+1)*k+j] += sum10;
                        (*res)[(i+2)*k+j] += sum20;
                }

                for(; j < k; j += 3) {
                        sum00 = sum01 = sum02 = sum10 = sum11 = sum12 = sum20 = sum21 = sum22 = 0;
		for(l = 0; l < m; l++) {
                                sum00 += a[i*m+l] * b[l*k+j];
                                sum01 += a[i*m+l] * b[l*k+j+1];
                                sum02 += a[i*m+l] * b[l*k+j+2];
                                sum10 += a[(i+1)*m+l] * b[l*k+j];
                                sum11 += a[(i+1)*m+l] * b[l*k+j+1];
                                sum12 += a[(i+1)*m+l] * b[l*k+j+2];
                                sum20 += a[(i+2)*m+l] * b[l*k+j];
                                sum21 += a[(i+2)*m+l] * b[l*k+j+1];
                                sum22 += a[(i+2)*m+l] * b[l*k+j+2];
                        }
                        (*res)[i*k+j] += sum00;
                        (*res)[i*k+j+1] += sum01;
                        (*res)[i*k+j+2] += sum02;
                        (*res)[(i+1)*k+j] += sum10;
                        (*res)[(i+1)*k+j+1] += sum11;
                        (*res)[(i+1)*k+j+2] += sum12;
                        (*res)[(i+2)*k+j] += sum20;
                        (*res)[(i+2)*k+j+1] += sum21;
                        (*res)[(i+2)*k+j+2] += sum22;
                }
        }
}

void subtract_vector(int n, const std::vector<double>& a,
    const std::vector<double>& b, std::vector<double>* c) {
    //std::vector<double> result(n);   
    for (int i = 0; i < n; ++i) {
        (*c)[i] = a[i] - b[i];
    }
}

double vector_norm(int n, const std::vector<double>& vec) {
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        norm += std::fabs(vec[i]);        
    }

    return norm;
}

void get_block(int i, int j, int n, int m, int k, int l,
 const std::vector<double>& matrix, std::vector<double>* block1) {
    int h = i < k ? m : l;
    int w = j < k ? m : l;
    
    int ind = 0;

    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            (*block1)[ind] = matrix[n * (m * i + p) + m * j + q];
            ind++;
        }
    }
}

void put_block(int i, int j, int n, int m, int k, int l, const std::vector<double>& block, std::vector<double>* matrix) {
    int h = i < k ? m : l;
    int w = j < k ? m : l;

    int ind = 0;
    for (int p = 0; p < h; ++p) {
        for (int q = 0; q < w; ++q) {
            (*matrix)[n * (m * i + p) + m * j + q] = block[ind];
            ind++;
        }
    }
}

void get_vector(int i, int m, int k, int l, const std::vector<double>& b, std::vector<double>* block3) {
    int length = i < k ? m : l;
    for (int p = 0; p < length; ++p) {
        (*block3)[p] = b[m*i + p]; 
    }
}

void put_vector(int i, int m, int k, int l, const std::vector<double>& b_i, std::vector<double>* b) {
    int length = i < k ? m : l;
    for (int p = 0; p < length; ++p) {
        (*b)[m*i + p] = b_i[p];
    }
}

// a -= b;
// n x m shape of both matrix
void subtract_matrix_inplace(int n, int m, std::vector<double>* a, const std::vector<double>& b) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            (*a)[m*i + j] -= b[m*i + j];           
        }
    }
}

void inverse_matrix(int m, std::vector<double>* matrix, std::vector<double>* identity,
 std::vector<int>* rows) {

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            (*identity)[i * m + j] = (i != j) ? 0 : 1;
        }
    }

    for (int k = 0; k < m; ++k) {
        (*rows)[k] = k;
    }

    for (int i = 0; i < m; ++i) {
        // выбираем максимальный элемент по столбцу.
        double max_elem = (*matrix)[(*rows)[i] * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs((*matrix)[(*rows)[j] * m + i]) > std::fabs(max_elem)) {
                max_elem = (*matrix)[(*rows)[j] * m + i];
                row_max_elem = j;
            }
        }

        std::swap((*rows)[i], (*rows)[row_max_elem]);

        // делим нашу строку на max_elem.
        double factor = 1 / max_elem;
        for (int s = 0; s < i; ++s) {
            (*identity)[(*rows)[i] * m + s] *= factor;
        }

        for (int s = i; s < m; ++s) {
            (*matrix)[(*rows)[i] * m + s] *= factor;
            (*identity)[(*rows)[i] * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -(*matrix)[(*rows)[k] * m + i];
            for (int p = 0; p < i + 1; ++p) {
                (*identity)[(*rows)[k] * m + p] += (*identity)[(*rows)[i] * m + p] * multiplier;
            }

            for (int p = i + 1; p < m; ++p) { 
                (*matrix)[(*rows)[k] * m + p] += (*matrix)[(*rows)[i] * m + p] * multiplier;
                (*identity)[(*rows)[k] * m + p] += (*identity)[(*rows)[i] * m + p] * multiplier;
            }
        }
    }

    for (int i = m - 1; i > 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            double multiplier = -(*matrix)[(*rows)[k] * m + i];
            for (int p = 0; p < m; ++p) { 
                (*identity)[(*rows)[k] * m + p] += (*identity)[(*rows)[i] * m + p] * multiplier;
            }
        }
    }

    // rows - настоящий порядок строк у обратной матрицы.
}

bool is_inv(int m, std::vector<double>* matrix, double a_norm) {
    std::vector<int> rows(m);

    for (int k = 0; k < m; ++k) {
        rows[k] = k;
    }

    for (int i = 0; i < m; ++i) {
        // выбираем максимальный элемент по столбцу.
        double max_elem = (*matrix)[rows[i] * m + i];
        int row_max_elem = i;
        for (int j = i + 1; j < m; ++j) {
            if (std::fabs((*matrix)[rows[j] * m + i]) - EPS * a_norm > std::fabs(max_elem)) {
                max_elem = (*matrix)[rows[j] * m + i];
                row_max_elem = j;
            }
        }

        std::swap(rows[i], rows[row_max_elem]);

        if (std::fabs(max_elem) < EPS * a_norm) {
            return false;    
        }

        // делим нашу строку на max_elem.
        double factor = 1 / max_elem;
        for (int s = i; s < m; ++s) {
            (*matrix)[rows[i] * m + s] *= factor;
        }

        for (int k = i + 1; k < m; ++k) {
            double multiplier = -(*matrix)[rows[k] * m + i];

            for (int p = i + 1; p < m; ++p) { 
                (*matrix)[rows[k] * m + p] += (*matrix)[rows[i] * m + p] * multiplier;
            }
        }
    }

    return true;
}
