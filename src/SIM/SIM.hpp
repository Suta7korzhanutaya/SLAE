#pragma once
#include <iostream>
#include <Matrix & Vector/CSR.hpp>
#include <Matrix & Vector/vector.hpp>
#include <vector>
std::vector<double> Simple(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tau, long double tolerance);

std::vector<double> Jacobi( CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tolerance);

std::vector<double> Gauss_Seidel(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tolerance);

std::vector<double> Chebyshev(const CSR& A, const std::vector<double> &b, const std::vector<double> &x0, long double tolerance, std::size_t res, double lambda_min, double lambda_max);

std::vector<double> Sym_Gauss_Seidel(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, double tolerance);

std::vector<double> Chebyshev_Sym_GS(const CSR& A,const std::vector<double>& b, const std::vector<double>& x0, double tolerance, double rho);

std::vector<double> iteration_Sym_GS(const CSR& A, const std::vector<double> &b, const std::vector<double> &y);
