#include <vector>
#include <cmath>
#include "inc.h"

double r1_eval(const std::vector<double>& matrix, const std::vector<double>& x,
    const std::vector<double>& b, std::vector<double>* c, std::vector<double>* d) {
    
    int n = x.size();
    
    matr_prod(n, n, 1, matrix, x, c);
    subtract_vector(n, *c, b, d);

    return vector_norm(*d) / vector_norm(b);
}

double r2_eval(const std::vector<double>& x) {
    int n = x.size();
    double r2 = 0;
    for (int i = 1; i <= n; ++i) {
        r2 += std::fabs(x[i - 1] - (i % 2));
    }

    return r2;
}
