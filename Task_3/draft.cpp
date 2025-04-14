#include <cmath>
#include "Dense.hpp"
#include "vector.hpp"

void Householder(const Dense& A) {
    std::size_t m = A.get_rows();
    std::size_t n = A.get_cols();
    Dense R = A;
    Dense Q = single(m);

    for (int h = 0; h < n; ++h) {
        std::vector<double> x(m - h);
        for (int i = h; i < m; ++i) {
            x[i - h] = R(i, h);
        }

        double norm_x = 0.0;
        for (int i = 0; i < m - h; ++i) {
            norm_x += x[i] * x[i];
        }
        
        norm_x = sign(x) * modul(x);

        std::vector<double> v(m - h);
        v[0] = x[0] + norm_x;
        for (int i = 1; i < m - h; ++i) {
            v[i] = x[i];
        }

        Dense H = single(m - h) - remain(v, v) * (2.0 / (norm_x * norm_x));
        Dense H_full = single(m);
        
        for (int i = h; i < m; ++i) {
            for (int j = h; j < m; ++j) {
                H_full(i, j) = H(i - h, j - h);
            }
        }

        R = H_full * R;
        Q = Q * H_full;
    }
}
