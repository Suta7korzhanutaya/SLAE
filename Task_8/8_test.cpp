#include <iostream>
#include <fstream>
#include <vector>

#include "../src/Gradient/GM.hpp"
#include "../src/GMRES/GMRES.hpp"
#include <chrono>

using high_resolution_clock = std::chrono::high_resolution_clock;

void get_time_Conjugate_Gradient(){
	std::ofstream file;
        file.open("res_CG.txt"); //файл для записи данных
	std::vector<double> a = {1250.5, 7.8, 5.2, 0.9, 3.1, 8.7, 568.4, 5.1, 12.4, 3.4, 4.7, 8.5, 1284.3, 25.4, 0.8, 8.6, 5.4, 2.3, 936.7, 6.2, 7.4, 1.5, 6.7, 4.5, 935.7};
	std::vector<double> x = {3.5, 1.7, 2.8, 12.5, 7.3};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		double alpha;
    		double beta;
		
	auto begin = high_resolution_clock::now(); // начало записи времени
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    
	    	std::vector<double> d = res;
    		std::vector<double> r = res; // начальная невязка
		
    		while(modul(res) > tolerance){
			alpha = (res * res) / (d * (A * d));
			y = y - alpha * d;
			res = A * y - b;
		    	  
		    	beta = (res * res) / (r * r);
			d = res + beta * d;
			r = res;  
		}	
	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	
	file.close();
}

void get_time_GMRES(){
	std::ofstream file;
        file.open("res_GMRES.txt"); //файл для записи данных
	std::vector<double> a = {1250.5, 7.8, 5.2, 0.9, 3.1, 8.7, 568.4, 5.1, 12.4, 3.4, 4.7, 8.5, 1284.3, 25.4, 0.8, 8.6, 5.4, 2.3, 936.7, 6.2, 7.4, 1.5, 6.7, 4.5, 935.7};
	std::vector<double> x = {3.5, 1.7, 2.8, 12.5, 7.3};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
	for(std::size_t m = 1; m < 100; m++){
		auto begin = high_resolution_clock::now(); // начало записи времени
		    	
		    	std::pair<std::vector<double>, double> result = Step(A, b, x0, m);
	    		while(result.second > tolerance){
	    			result = Step(A, b, result.first, m);
	    		}
	    			
		auto end = high_resolution_clock::now(); // конец записи времени
		auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
		file << time.count() << std::endl;
	}
	file.close();
}

int main(){
    get_time_Conjugate_Gradient();
    get_time_GMRES();  
}
