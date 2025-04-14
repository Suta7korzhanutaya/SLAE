#pragma once
#include <iostream>
#include <Matrix & Vector/CSR.hpp>
#include <Matrix & Vector/vector.hpp>
#include <vector>

void Iteration(const CSR& A, Dense& K, Dense& H, const std::size_t i);
void Givens_rot(Dense& H, std::vector<std::pair<double, double>>& rot, const std::size_t i);

std::vector<double> Back_Gauss(const Dense& A, std::vector<double> b) ;

std::pair<std::vector<double>, double> Step(const CSR& A, const std::vector<double>& b, std::vector<double> x0, const std::size_t m);


inline std::vector<double> GMRES(CSR& A, const std::vector<double>& b,  const std::vector<double>& x0, const std::size_t m, long double tolerance){
    std::pair<std::vector<double>, double> result = Step(A, b, x0, m);
    
    while(result.second > tolerance){
    	result = Step(A, b, result.first, m);
    }
    
    return result.first;
}
