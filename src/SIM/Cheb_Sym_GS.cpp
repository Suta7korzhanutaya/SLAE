#include "SIM.hpp"

std::vector<double> Chebyshev_Sym_GS(const CSR& A,const std::vector<double>& b, const std::vector<double>& x0, double tolerance, double rho) {
            std::vector<double> y_m = x0;
            std::vector<double> y_p = iteration_Sym_GS(A, b, y_m);
            
            for (std::size_t i = 0; i < y_p.size(); i++)
            {
            	std::cout<<y_p[i]<<std::endl;
       
            }
            
            std::cout<< "nmsjdvn" <<std::endl;
            
            std::vector<double> res = A * y_p - b;

            double mu_m = 1;
            double mu = 1 / rho;
            double mu_p;

            while (modul(res) > tolerance) {
                mu_p = 2 * mu / rho - mu_m;

                y_p = iteration_Sym_GS(A, b, y_p);

	

                for (std::size_t i = 0; i < y_p.size(); i++) {
                    y_p[i] = 2 * mu * y_p[i] / ( rho * mu_p) - mu_m * y_m[i] / mu_p;
                }
                
                for (std::size_t i = 0; i < y_p.size(); i++)
            {
            	std::cout<<y_p[i]<<std::endl;
            }

                y_m = y_p;
                mu_m = mu;
                mu = mu_p;

		std::cout<< mu_m <<std::endl;
		std::cout<< mu <<std::endl;
	

                res = A * y_p - b;
            }
            return y_p;
        }



std::vector<double> iteration_Sym_GS(const CSR& A, const std::vector<double> &b, const std::vector<double> &y) {
            std::vector<double> x = y;
            double diag;

            for (std::size_t i = 0; i < x.size(); ++i) {
                x[i] = b[i];
                for (std::size_t k = A.get_rows()[i]; k < A.get_rows()[i+1]; k++) {
                    if (A.get_cols()[k] == i) {
                        diag = A(i,i);
                        continue;
                    }
                    x[i] -= A.get_values()[k] * x[A.get_cols()[k]];
                }
                x[i] /= diag;
            }

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
            return x;
}
