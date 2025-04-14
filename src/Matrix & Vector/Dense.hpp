#pragma once

#include <iostream>
#include <vector>
#include "vector.hpp"

class Dense 
{
	private:
	    std::size_t rows, cols;
	    std::vector<double> data;
	public:
	    // конструктор класса 
	    Dense(const std::vector<double>& data, std::size_t rows, std::size_t cols);
	    
	    // доступ к элементу матрицы по индексам для чтения
	    double operator()(std::size_t i, std::size_t j) const;
	    
	    // количество столбцов
	    std::size_t get_cols() const;
	    
	    // количество строк
	    std::size_t get_rows() const;
	    
	    // сложение двух матриц
	    Dense operator+(const Dense& other) const;
	    
	    // умножение матрицы на число
	    Dense operator*(double alpha) const;
	    
	    // умножение матрицы на вектор
	    std::vector<double> operator*(const std::vector<double>& vec) const;
	    
	    // умножение двух матриц
	    Dense operator*(const Dense& other) const;
	    
	    // вывод матрицы
	    void print() const;
	    
	    Dense(const double elem, const std::size_t rows, const std::size_t cols);
	    
	    double& get_element(std::size_t i,  std::size_t j);
	    void change_column(const std::vector<double>& vec, std::size_t i);
	    
	    const std::vector<double> get_column(std::size_t num) const;
	    const std::vector<double> get_row(std::size_t num) const;
	    std::vector<double> get_elements(const std::vector<double>& data) const;
	    void single(std::size_t n);
	    void swap_elements(std::size_t i, std::size_t j, double num);
	    Dense transpon() const;
	    
	    bool operator == (const Dense & other);
};
