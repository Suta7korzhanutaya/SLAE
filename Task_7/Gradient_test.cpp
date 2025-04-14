#include <iostream>
#include <fstream>
#include <vector>

#include "../src/SIM/SIM.hpp"
#include "../src/Gradient/GM.hpp"
#include <chrono>

using high_resolution_clock = std::chrono::high_resolution_clock;

void get_time_Simple_iteration(){
	std::ofstream file;
        file.open("res_Simple.csv"); //файл для записи данных
	std::vector<double> a = {1250.5, 7.8, 5.2, 0.9, 3.1, 8.7, 568.4, 5.1, 12.4, 3.4, 4.7, 8.5, 1284.3, 25.4, 0.8, 8.6, 5.4, 2.3, 936.7, 6.2, 7.4, 1.5, 6.7, 4.5, 935.7};
	std::vector<double> x = {3.5, 1.7, 2.8, 12.5, 7.3};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		double tau = 0.00108;
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
		
		std::size_t n = 0;
    		auto begin = high_resolution_clock::now(); // начало записи времени
    		
    		while(modul(res) > tolerance){
			y = y - tau * res;
	    		res = A * y - b;
		        
		        auto now = high_resolution_clock::now();
		    	auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    	file << modul(res) << " " << n << " "<< t.count() << std::endl;
		    	n++;
		        
         	}
		auto end = high_resolution_clock::now(); // конец записи времени
		auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
		file << time.count() << std::endl;
	file.close();
}



void get_time_Chebyshev(){
	std::ofstream file;
        file.open("res_Chebyshev.csv"); //файл для записи данных
	std::vector<double> a = {1250.5, 7.8, 5.2, 0.9, 3.1, 8.7, 568.4, 5.1, 12.4, 3.4, 4.7, 8.5, 1284.3, 25.4, 0.8, 8.6, 5.4, 2.3, 936.7, 6.2, 7.4, 1.5, 6.7, 4.5, 935.7};
	std::vector<double> x = {3.5, 1.7, 2.8, 12.5, 7.3};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	std::vector<double> b = A * x;  // вектор свободных членов
	double tolerance = pow(10,-12);
	double lambda_max = 1285.395;
	double lambda_min = 568.055;
	std::size_t res = 5;
	
		double R = static_cast<double>(res);
		double n = pow(2,R); //количество прогонок
		
		std::size_t N = static_cast<std::size_t>(n);
		
		std::vector<std::size_t > res_num(N,0);
		std::size_t P = static_cast<std::size_t>(pow(2,R-1));
	    	res_num[P] = 1;
	    	
	    	for(std::size_t  i = 2; i <= res; i++){
	    		int I = static_cast<int>(i);
	    		std::size_t P1 = static_cast<std::size_t>(pow(2,R-I));
			std::size_t k = P1;
			for(std::size_t  j = 0; j < N; j+= 2*k){
				std::size_t P2 = static_cast<std::size_t>(pow(2,I));
		    		res_num[j+k] = P2 - res_num[j] - 1;
			}
	    	}
	    
	    	double cos_a = cos(M_PI/n);
	    	double sin_a = sin(M_PI/n);
	    	double sin_b = sin(M_PI/(2*n));
	    
	    	std::vector<double> z(N,0);
	    	z[0] = sin(M_PI/(2*n));

	    	for(std::size_t  i = 1; i < N; i++){
			z[i] = z[i-1]*cos_a - sin_b*sin_a;
			sin_b = z[i-1]*cos_a + sin_b*sin_a;
	       	 	z[i-1] = z[i-1]*(lambda_max-lambda_min)/2 + (lambda_max+lambda_min)/2;
	    	}
	    
	    	z[N-1] = z[N-1]*(lambda_max-lambda_min)/2 + (lambda_max+lambda_min)/2;



    	std::vector<double> y = x0; // x0 - начальное приближение
    	std::vector<double> r = A * y - b;
   
   	std::size_t T = 0;
    	auto begin = high_resolution_clock::now(); // начало записи времени
   
    	while(modul(r) > tolerance){
        	for(std::size_t i = 0; i < res_num.size(); i++) {
            		y = y - r * (1 / z[res_num[i]]);
            		r = A * y - b;
            		auto now = high_resolution_clock::now();
		    	auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    	file << modul(r) << " " << T << " "<< t.count() << std::endl;
		    	T++;
        	}
    	}
    	
    	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	file.close();
}

void get_time_Gradient(){
	std::ofstream file;
        file.open("res_Gradient.csv"); //файл для записи данных
	std::vector<double> a = {1250.5, 7.8, 5.2, 0.9, 3.1, 8.7, 568.4, 5.1, 12.4, 3.4, 4.7, 8.5, 1284.3, 25.4, 0.8, 8.6, 5.4, 2.3, 936.7, 6.2, 7.4, 1.5, 6.7, 4.5, 935.7};
	std::vector<double> x = {3.5, 1.7, 2.8, 12.5, 7.3};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	double tau; 
		
		std::size_t n = 0;
		auto begin = high_resolution_clock::now(); // начало записи времени
    		while(modul(res) > tolerance){
			tau = (res * res) / ( res * ( A * res) );
		    	y = y - tau * res;
		    	res = A * y - b;
		    
		    auto now = high_resolution_clock::now();
		    auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    file << modul(res) << " " << n << " "<< t.count() << std::endl;
		    	n++;
		}	
	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	
	file.close();
}

void get_time_Conjugate_Gradient(){
	std::ofstream file;
        file.open("res_Con_Gradient.csv"); //файл для записи данных
	std::vector<double> a = {1250.5, 7.8, 5.2, 0.9, 3.1, 8.7, 568.4, 5.1, 12.4, 3.4, 4.7, 8.5, 1284.3, 25.4, 0.8, 8.6, 5.4, 2.3, 936.7, 6.2, 7.4, 1.5, 6.7, 4.5, 935.7};
	std::vector<double> x = {3.5, 1.7, 2.8, 12.5, 7.3};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 5, 5);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	
	    	double alpha;
    		double beta;
	    	
	    	std::vector<double> d = res;
    		std::vector<double> r = res; // начальная невязка
		
		std::size_t n = 0;
		auto begin = high_resolution_clock::now(); // начало записи времени
    		while(modul(res) > tolerance){
			alpha = (res * res) / (d * (A * d));
			y = y - alpha * d;
			res = A * y - b;
		    	  
		    	beta = (res * res) / (r * r);
			d = res + beta * d;
			r = res;  	
		    
		    auto now = high_resolution_clock::now();
		    auto t = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
		    file << modul(res) << " " << n << " "<< t.count() << std::endl;
		    	n++;
		}	
	auto end = high_resolution_clock::now(); // конец записи времени
	auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
	file << time.count() << std::endl;
	
	file.close();
}


int main(){
    get_time_Simple_iteration();
    get_time_Chebyshev();
    get_time_Gradient();
    get_time_Conjugate_Gradient();  
}
