#pragma once

#include <iostream>
#include <vector>
#include "Dense.hpp"
#include "vector.hpp"


class CSR {
private:
    std::vector<double> values;
    std::vector<unsigned int> cols;
    std::vector<unsigned int> rows;

public:  
    // конструктор класса 
    CSR(const std::vector<double>& values, const std::vector<unsigned int>& cols, const std::vector<unsigned int>& rows);
    
    // конструктор класса 
    CSR(const std::vector<double>& data, std::size_t r, std::size_t c);
    
    // доступ к элементу матрицы для чтения
    double operator()(std::size_t i, std::size_t j) const;
    
    // умножение матрицы на число
    CSR operator*(double scalar) const;
    
    // умножение матрицы на вектор
    std::vector<double> operator*(const std::vector<double>& vec) const;
    
    // вывод матрицы в нужном формате
    void print() const;
    
    std::vector<double> Jakobi_reverse_diag() const;
    std::vector<double> Jakobi_non_diagonal(const std::vector<double>& vec) const;

    // получение векторов
    const std::vector<double>& get_values() const {
    	return values;
    }
    
    const std::vector<unsigned int>& get_cols() const {
    	return cols;
    }
    
    const std::vector<unsigned int>& get_rows() const {
    	return rows;
    }
};
