#include "SIM.hpp"
#include <cmath>

std::vector<double> Chebyshev(const CSR& A, const std::vector<double> &b, const std::vector<double> &x0, long double tolerance, std::size_t res, double lambda_min, double lambda_max){
                             
    std::size_t n = pow(2,res); //количество прогонок
    std::vector<std::size_t > res_n(n,0);
    res_n[pow(2,res-1)] = 1;
    for(std::size_t  i = 2; i <= res; i++){
        std::size_t k = pow(2,res-i);
        for(std::size_t  j = 0; j < n; j+= 2*k){
            res_n[j+k] = pow(2,i) - res_n[j] - 1;
        }
    }
    
    double cos_a = cos(M_PI/n);
    double sin_a = sin(M_PI/n);
    double sin_b = sin(M_PI/(2*n));
    
    std::vector<double> x(n,0);
    x[0] = sin(M_PI/(2*n));

    for(std::size_t  i = 1; i < n; i++){
        x[i] = x[i-1]*cos_a - sin_b*sin_a;
        sin_b = x[i-1]*cos_a + sin_b*sin_a;
        x[i-1] = x[i-1]*(lambda_max-lambda_min)/2 + (lambda_max+lambda_min)/2;
    }
    
    x[n-1] = x[n-1]*(lambda_max-lambda_min)/2 + (lambda_max+lambda_min)/2;

    std::vector<double> y = x0;
    std::vector<double> r = A * y - b;
    
    while(modul(r) > tolerance){
        for(std::size_t i = 0; i < res_n.size(); i++) {
            y = y - r * (1 / x[res_n[i]]);
            r = A * y - b;
        }
    }
    return y;
}
