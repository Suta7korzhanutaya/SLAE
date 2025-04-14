#include <iostream>
#include <fstream>
#include <vector>

#include <random>
#include <chrono>

#include "../src/Matrix & Vector/Dense.hpp"
#include "../src/Matrix & Vector/CSR.hpp"

using high_resolution_clock = std::chrono::high_resolution_clock;

void get_time_CSR(){
	std::srand(0);  //генерация случайных чисел
	std::vector<double> a; //вектор для генерации матрицы
    	std::vector<double> b; //вектор для умножения на матрицу
	double val; //значение из вектора a
    std::ofstream file;
    file.open ("time_CSR.csv"); //файл для записи данных
	for(size_t i = 10; i < 2110; i += 10){
		for(size_t j = 0; j < i * i; j++){
			val = ((std::rand() % 5 == 0) ? std::rand() : 0); 
			a.push_back(val);
		}
		Dense D(a, i, i);
		CSR A(D);

		for(size_t j = 0; j <= i; j++){
			b.push_back(std::rand());
		}
        
        file << i << ", ";
		auto begin = high_resolution_clock::now();
		std::vector<double> result = A * b;
		auto end = high_resolution_clock::now();
		auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
        file << time.count() << std::endl; 
    }
    file.close();
}



void get_time_Dense(){
	std::srand(0); //генерация случайных чисел
	std::vector<double> a; //вектор для генерации матрицы
    	std::vector<double> b; //вектор для умножения на матрицу
	double val; //значение из вектора a
    std::ofstream file;
    file.open ("time_Dense.csv"); //файл для записи данных
	for(size_t i = 10; i < 2110; i += 10){
		for(size_t j = 0; j < i * i; j++){
			val = ((std::rand() % 5 == 0) ? std::rand() : 0); 
			a.push_back(val);
		}
		Dense A(a, i, i);

		for(size_t j = 0; j <= i; j++){
			b.push_back(std::rand());
		}
        
        file << i << ", ";
		auto begin = high_resolution_clock::now();
		std::vector<double> result = A * b;
		auto end = high_resolution_clock::now();
		auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
        file << time.count() << std::endl; 
    }
    file.close();
}


int main(){
    get_time_Dense();
    get_time_CSR();
}
