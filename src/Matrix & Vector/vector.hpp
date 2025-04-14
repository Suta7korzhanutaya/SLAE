#pragma once
#include <vector>
#include <iostream>
#include <cmath>

// Умножение векторв на вектор
double operator*(const std::vector<double>& a, const std::vector<double>& b);

// Умножение вектора на число
std::vector<double> operator*(const std::vector<double>& a, const double& alpha);

// Сложение векторов
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

// Умножение числа на вектор
std::vector<double> operator*(const double& alpha, const std::vector<double>& a);

// Вычитание векторов
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

void print(const std::vector<double>& result) ;

std::vector<double> operator/(const std::vector<double>& a, const double& alpha);

std::vector<double> transposed(const std::vector<double>& x);

double modul(const std::vector<double>& x);
double sq(const std::vector<double>& x);

int sign(const std::vector<double>& x);
std::vector<double> mul(const std::vector<double>& a, const std::vector<double>& b);
