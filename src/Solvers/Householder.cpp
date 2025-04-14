#include <utility>
#include <cmath>
#include <Matrix & Vector/Dense.hpp>
#include <Matrix & Vector/vector.hpp>

std::pair<Dense, Dense> Householder(const Dense& A) {
    Dense R = A;
    std::size_t m = R.get_rows();
    std::size_t b = R.get_cols();
    std::vector<double> e(m, 0); // базисный вектор e_1
    e[0] = 1;
    std::vector<double> v = R.get_column(0) + sign(R.get_column(0)) * modul(R.get_column(0)) * e;
    
    std::vector<double> I(b * m, 0);
    	for (std::size_t k = 0; k < b; ++k) {
        	I[k * b + k] = 1.0;
        }
    
    Dense Q(I, m, b);

    for (unsigned int n = 0; n < R.get_cols(); n++) {
        std::size_t size = m - n;
        std::vector<double> x = R.get_column(n);
        std::vector<double> vec = x + sign(x) * modul(x) * e;

        std::vector<double> prod = x - 2 * vec * ((vec * x) / std::pow(modul(vec), 2));
        R.swap_elements(n, n, prod[0]);

        for (unsigned int  i = 1 + n; i < b; i++) {
            R.swap_elements(i, n, 0);
        }

        for (unsigned int  j = 1 + n; j < b; j++) {
            std::vector<double> a = R.get_column(j);
            x = std::vector<double>(a.begin() + n, a.end());
            prod = x - 2 * vec * ((vec * x) / std::pow(modul(vec), 2));

            for (unsigned int i = 0; i < size; i++) {
                R.swap_elements(i + n, j, prod[i]);
            }
        }

        for (unsigned int i = 0; i < b; i++) {
            std::vector<double> a = Q.get_row(i);
            x = std::vector<double>(a.begin() + n, a.end());

            for (unsigned int j = 0; j < size; j++) {
                prod = x - 2 * vec * ((vec * x) / std::pow(modul(vec), 2));
                Q.swap_elements(i, j + n, prod[j]);
            }
        }
    }

    return std::make_pair(Q, R);
}
