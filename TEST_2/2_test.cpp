#include <iostream>
#include <fstream>
#include <vector>

#include "../src/SIM/SIM.hpp"
#include "../src/Gradient/GM.hpp"
#include <chrono>

void get_time_1(){
	std::ofstream file;
        file.open("res_1.csv"); //файл для записи данных
	std::vector<double> a = {6, 0, 0, 0, 0, 7, 0, 0, 0, 0, 8, 0, 0, 0, 0, 9};
	std::vector<double> x = {1, 0.857, 0.75, 0.67};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 4, 4);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		double tau = 0.9 * 2/9;
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
		
    		while(modul(res) > tolerance){
			y = y - tau * res;
	    		res = A * y - b;
		        
		    	file << y[0] <<" "<< y[3] << std::endl;   
         	}
	file.close();
}


void get_time_2(){
	std::ofstream file;
        file.open("res_2.csv"); //файл для записи данных
	std::vector<double> a = {6, 0, 0, 0, 0, 7, 0, 0, 0, 0, 8, 0, 0, 0, 0, 9};
	std::vector<double> x = {1, 0.857, 0.75, 0.67};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 4, 4);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		double tau = 0.133;
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
		
    		while(modul(res) > tolerance){
			y = y - tau * res;
	    		res = A * y - b;
		        
		        
		    	file << y[0] <<" "<< y[3] << std::endl;   
         	}
	file.close();
}




void get_time_Chebyshev(){
	std::ofstream file;
        file.open("res_4.csv"); //файл для записи данных
	std::vector<double> a = {6, 0, 0, 0, 0, 7, 0, 0, 0, 0, 8, 0, 0, 0, 0, 9};
	std::vector<double> x = {1, 0.857, 0.75, 0.67};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 4, 4);

	std::vector<double> b = A * x;  // вектор свободных членов
	double tolerance = pow(10,-12);
	double lambda_max = 9;
	double lambda_min = 6;
	std::size_t res = 8;
	
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
   
    
    	while(modul(r) > tolerance){
        	for(std::size_t i = 0; i < res_num.size(); i++) {
            		y = y - r * (1 / z[res_num[i]]);
            		r = A * y - b;
            		file << y[0] <<" "<< y[3] << std::endl;   
        	}
    	}
    	
	file.close();
}

void get_time_Gradient(){
	std::ofstream file;
        file.open("res_3.csv"); //файл для записи данных
	std::vector<double> a = {6, 0, 0, 0, 0, 7, 0, 0, 0, 0, 8, 0, 0, 0, 0, 9};
	std::vector<double> x = {1, 0.857, 0.75, 0.67};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 4, 4);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	double tau; 
		
    		while(modul(res) > tolerance){
			tau = (res * res) / ( res * ( A * res) );
		    	y = y - tau * res;
		    	res = A * y - b;
		    
		    file << y[0] <<" "<< y[3] << std::endl;   
		}	

	file.close();
}

void get_time_Conjugate_Gradient(){
	std::ofstream file;
        file.open("res_5.csv"); //файл для записи данных
	std::vector<double> a = {6, 0, 0, 0, 0, 7, 0, 0, 0, 0, 8, 0, 0, 0, 0, 9};
	std::vector<double> x = {1, 0.857, 0.75, 0.67};
	std::vector<double> x0 = {0, 0, 0, 0, 0};
	
	CSR A(a, 4, 4);
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-12);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	
	    	double alpha;
    		double beta;
	    	
	    	std::vector<double> d = res;
    		std::vector<double> r = res; // начальная невязка
		
		
    		while(modul(res) > tolerance){
			alpha = (res * res) / (d * (A * d));
			y = y - alpha * d;
			res = A * y - b;
		    	  
		    	beta = (res * res) / (r * r);
			d = res + beta * d;
			r = res;  	
		   
		   file << y[0] <<" "<< y[3] << std::endl;   
		}	
	
	file.close();
}


int main(){
    get_time_1();
    get_time_2();
    get_time_Gradient();
    get_time_Chebyshev();
    get_time_Conjugate_Gradient();  
}
