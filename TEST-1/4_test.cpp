#include <iostream>
#include <fstream>
#include <vector>

#include "../src/SIM/SIM.hpp"
#include <chrono>


void get_time_Simple_iteration(){
	std::ofstream file;
        file.open("res.csv"); //файл для записи данных
	std::vector<double> a = {56, 0, 5, 4, 8, 76, 7, 9, 0, 1, 56, 0, 0, 8, 6, 79};
	std::vector<double> x0 = {0, 0, 0, 0};
	
	CSR A(a, 4, 4);
	for(double tau = 0.001; tau <= 0.023; tau += 0.0001){ 
		std::vector<double> b = {1, 1, 1, 1};  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
		
		std::size_t n = 0;
    		while(modul(res) > tolerance){
			y = y - tau * res;
	    		res = A * y - b;
	    		n++;
         	}
         	file << tau << " " << n << std::endl;
	}    	
	file.close();
}

void get_time_Jacobi(){
	std::ofstream file;
        file.open("res1.csv"); //файл для записи данных
	std::vector<double> a = {97, 0, 4, 0, 0, 60, 7, 0, 0, 0, 69, 0, 2, 3, 5, 50};
	std::vector<double> x0 = {0, 0, 0, 0};
	
	CSR A(a, 4, 4);

	std::vector<double> b = {5, 7, 7, 8};  // вектор свободных членов
		
	std::vector<double> y = x0; // x0 - начальное приближение
	std::vector<double> res = A * y - b;
	std::vector<double> LU;
	std::vector<double> t = res;		
	std::size_t n = 0;
	
    	while (modul(res) > (modul(t))/1000){
		LU = A.Jakobi_non_diagonal(y);
		y = mul(A.Jakobi_reverse_diag(), (b - LU));
		res = A * y - b;
		    	
		file << modul(res) << " " << n << std::endl; 	
		n++;
	} 	

	file.close();
}

int main(){
    get_time_Jacobi();
    get_time_Simple_iteration();
}
