#include "GMRES.hpp"

std::vector<double> Back_Gauss(const Dense& A, std::vector<double> b) {
    std::vector<double> result = b;
    std::size_t n = b.size() - 1;
    result[n] = b[n] / A(n, n);
    
    for (std::size_t i = n; i > 0; i--) {
        double norm = b[i - 1];
        
        for (std::size_t k = n; k > i - 1; k--) {
        	norm -= A(i - 1, k) * result[k];
        }
        
        result[i - 1] = norm / A(i - 1, i - 1);
    }
    
    return result;
}
