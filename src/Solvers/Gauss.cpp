#include <Matrix & Vector/Dense.hpp>
#include "Householder.cpp"

std::vector<double> Gauss(const Dense& A, const std::vector<double>& b) {
	Dense Q = Householder(A).first;
    	Dense R = Householder(A).second;
    	std::vector<double> d = Q.transpon() * b;
    	
    	std::size_t n = d.size();
    	std::vector<double> x(n, 0);
    	
    	int t = static_cast<int>(d.size()); //интовый итератор для цикла	
	double rev;
	x[n] = R(n, n);
    	for (int i = t - 1; i >= 0; i--) {
            rev = 0;
            std::size_t k = static_cast<std::size_t>(i); //size_t итератор для матрицы
            for (int j = t - 1; j > i; j--) {
                std::size_t h = static_cast<std::size_t>(j); //size_t итератор для матрицы
                rev += R(k, h) * x[h];
            }
            x[k] = (d[k] - rev) / R(k, k); 
        }
        return x;
}
