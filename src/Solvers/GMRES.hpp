#pragma once
#include <iostream>
#include <Matrix & Vector/CSR.hpp>
#include <Matrix & Vector/vector.hpp>
#include <Matrix & Vector/Dense.hpp>
#include <vector>

std::vector<double> Back_Gauss(const Dense& A, std::vector<double> b);
void Iteration(const CSR& A, Dense& K, Dense& H, const std::size_t i);
void Givens_rot(Dense& H, std::vector<std::pair<double, double>>& rot, const std::size_t i);

std::pair<std::vector<double>, double> Step(const CSR& A, const std::vector<double>& b, std::vector<double> x0, const std::size_t m);

/*
std::vector<double> Back_Gauss(const Dense& A, std::vector<double> b) {
    std::vector<double> result = b;
    std::size_t n = b.size() - 1;
    result[n] = b[n] / A(n, n);
    
    for (std::size_t i = n; i > 0; i--) {
        double norm = b[i - 1];
        
        for (std::size_t k = n; k > i - 1; k--) {
        	norm -= A(i - 1, k) * result[k];
        }
        
        result[i - 1] = norm / A(i - 1, i - 1);
    }
    
    return result;
}



void Iteration(const CSR& A, Dense& K, Dense& H, const std::size_t i){
    std::vector<double> v = A * K.get_column(i);
    double h;
    
    //заполняем матрицу Хессинберга H
    for (std::size_t k = 0; k <= i; k++){
        h = v * K.get_column(k);
        v = v - h * K.get_column(k);
        H.get_element(k, i) = h;
    }

    h = modul(v);
    H.get_element(i + 1, i) = h;
    v = 1 / h * v;
    K.change_column(v, i + 1);
}



void Givens_rot(Dense& H, std::vector<std::pair<double, double>>& rot, const std::size_t i){
    //применяем все предыдущие вращения к новому столбцу H
    for (std::size_t k = 0; k < i; k++){
        double h = rot[k].first * H(k, i) - rot[k].second * H(k + 1, i);
        double h_n = rot[k].first * H(k+1, i) + rot[k].second * H(k, i);
        H.get_element(k, i) = h;
        H.get_element(k + 1, i) = h_n;
    }
    
    //новое вращение к столбцу H
    double norm = std::pow(H(i, i)* H(i, i) + H(i + 1, i) * H(i+1, i) , 0.5);
    double cos_th = H(i, i) / norm;
    double sin_th = -H(i + 1, i) / norm;
  
    double h  = cos_th * H(i, i) - sin_th * H(i + 1, i);
    double h_n = sin_th * H(i, i) + cos_th * H(i + 1, i);

    H.get_element(i, i) = h;
    H.get_element(i + 1, i) = h_n;

    rot[i] = (std::make_pair(cos_th, sin_th));
}



std::pair<std::vector<double>, double> Step(const CSR& A, const std::vector<double>& b, std::vector<double> x0, const std::size_t m){

    Dense H(0, m + 1, m); //матрица Хессинберга
    Dense K(0, x0.size(), m + 1); //базис Крылова

    std::vector<std::pair<double, double>> rot(m);

    std::vector<double> r0 = A * x0 - b; //начальная невязка
    std::vector<double> e(m + 1, 0);
    e[0] = modul(r0);

    K.change_column(1/modul(r0) * r0, 0); //первый вектор Крылова - нормированная невязка

    for(std::size_t i = 0; i < m ; i++){
        Iteration(A, K, H, i);
        Givens_rot(H, rot, i);
       
        double e_m = rot[i].first * e[i] - rot[i].second * e[i+1];
        double e_p = rot[i].second * e[i] + rot[i].first * e[i+1];

        e[i] = e_m;
        e[i + 1] = e_p;   
    }

    std::vector<double> vec(e.begin(), e.begin() + m);
    
    std::vector<double> y = Back_Gauss(H, vec);
    std::vector<double> result(x0.size(), 0);
    
    for (std::size_t k = 0; k < y.size(); k++) {
    	result = result + y[k] * K.get_column(k); //ЛК базиса Крылова
    }
    
    double ans = (e[m] < 0) ? -e[m] : e[m]; //абсолютное значение элемента
    return {x0 - result, ans};
}
*/

std::vector<double> GMRES(const CSR& A, const std::vector<double>& b,  const std::vector<double>& x0, const std::size_t m, long double tolerance){
    std::pair<std::vector<double>, double> result = Step(A, b, x0, m);
    
    while(result.second > tolerance) {
    	result = Step(A, b, result.first, m);
    }
    
    return result.first;
}
