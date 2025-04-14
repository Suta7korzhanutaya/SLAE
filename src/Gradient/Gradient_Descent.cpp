#include "GM.hpp"

std::vector<double> Gradient(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, long double tolerance) {
    	std::vector<double> x = x0; // x0 - начальное приближение
    	std::vector<double> res = A * x - b;
    	double tau;
    
    	while (modul(res) > tolerance){
    		tau = (res * res) / ( res * ( A * res) );
	    	x = x - tau * res;
	    	res = A * x - b;
    	}   	
        return x;
}
