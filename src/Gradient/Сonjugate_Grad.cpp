#include "GM.hpp"

std::vector<double> Conjugate_Gradient(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, long double tolerance) {
    	std::vector<double> x = x0; // x0 - начальное приближение
    	std::vector<double> res = A * x - b;
    	std::vector<double> d = res;
    	std::vector<double> r = res; // начальная невязка
    	
    	double alpha;
    	double beta;

    	while (modul(res) > tolerance){
	    alpha = (res * res) / (d * (A * d));
	    x = x - alpha * d;
	    res = A * x - b;
    	  
    	    beta = (res * res) / (r * r);
	    d = res + beta * d;
	    r = res;  	
    	}
    	
        return x;
}

