#pragma once
#include <iostream>
#include <Matrix & Vector/CSR.hpp>
#include <Matrix & Vector/vector.hpp>
#include <vector>

std::vector<double> Gradient(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, long double tolerance);

std::vector<double> Conjugate_Gradient(const CSR& A, const std::vector<double>& b, const std::vector<double>& x0, long double tolerance);
