#include <iostream>
#include <fstream>
#include <vector>

#include "../src/SIM/SIM.hpp"
#include <chrono>


void get_time_Simple_iteration(){
	std::ofstream file;
        file.open("res_Simple.csv"); //файл для записи данных
	std::vector<double> a = {56, 0, 5, 4, 8, 76, 7, 9, 0, 1, 56, 0, 0, 8, 6, 79};
	std::vector<double> x0 = {0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = {1, 1, 1, 1};  // вектор свободных членов
		double tolerance = pow(10,-12);
		double tau = 71,4745;
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
		
		std::size_t n = 0;
    		while(modul(res) > tolerance){
			y = y - tau * res;
	    		res = A * y - b;
		      
		    	file << modul(res) << " " << n << std::endl;
		    	n++;
         	}
	file.close();
}

int main(){
    get_time_Simple_iteration();
}
