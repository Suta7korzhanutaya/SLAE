#include "SIM.hpp"

std::vector<double> Simple(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tau, long double tolerance) {
    	std::vector<double> x = x0; // x0 - начальное приближение
    	std::vector<double> res = A * x - b;
    
    	while (modul(res) > tolerance){
	    	x = x - tau * res;
	    	res = A * x - b;
    	}   	
        return x;
}
