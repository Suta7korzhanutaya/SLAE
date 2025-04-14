#include "SIM.hpp"

std::vector<double> Sym_Gauss_Seidel(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tolerance) {
    	std::vector<double> x = x0; // x0 - начальное приближение
    	std::vector<double> res = A * x - b;
    	double diag; // диагональный элемент
    	
    	while (modul(res) > tolerance){
		for (std::size_t i = 0; i < x.size(); ++i) {
		    x[i] = b[i];
		    for(std::size_t k = A.get_rows()[i]; k < A.get_rows()[i+1]; k++){
		       if (A.get_cols()[k] == i){
		            diag = A(i,i);
		            continue;
		       }
		       x[i] -= A.get_values()[k] * x[A.get_cols()[k]];
		    }
		    x[i] /= diag;
		}
		
		for (unsigned i = 0; i < x.size(); i++) {
		    std::cout << x[i] << " ";
		}
		std::cout << std::endl;
		
		for (std::size_t i = x.size(); i > 0;) {
		    --i;
		    x[i] = b[i];
		    for(std::size_t k = A.get_rows()[i]; k < A.get_rows()[i+1]; k++){
		       if (A.get_cols()[k] == i){
		            diag = A(i,i);
		            continue;
		       }
		       x[i] -= A.get_values()[k] * x[A.get_cols()[k]];
		    }
		    x[i] /= diag;
		}
		
		for (unsigned i = 0; i < x.size(); i++) {
		    std::cout << x[i] << " ";
		}
		std::cout << std::endl;
		
	    	res = A * x - b;
    	}	
        return x;
}
