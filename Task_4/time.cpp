#include <iostream>
#include <fstream>
#include <vector>

#include <random>
#include <chrono>

#include "../src/SIM/Jacobi.cpp"
#include "../src/SIM/Gauss-Seidel.cpp"
#include "../src/SIM/Simple_iteration.cpp"

using high_resolution_clock = std::chrono::high_resolution_clock;

void get_time_Jacobi(){
	std::srand(0);  //генерация случайных чисел
	double val; //значение из вектора a
	std::ofstream file;
        file.open("res_Jacobi.txt"); //файл для записи данных
	std::vector<double> a(10000); //вектор для генерации матрицы
	std::vector<double> x(100);
	std::vector<double> x0(100);
	
		for(size_t i = 0; i < 100; i++){
		   for(size_t j = 0; j < 100; j++){
			val = (std::rand()% (5000 - 10 + 1)) + 10; 
			if (i != j) a[i*100 + j] = val;
		   }
		}
		
		for(size_t i = 0; i < 100; i++){
		   double sum = 0.0;
		   for(size_t j = 0; j < 100; j++){
			if (i != j){
				sum += a[i*100 + j];
			}
		   }
		   a[i*100 + i] = sum + (std::rand()% (5000 - 10 + 1)) + 10;
		}
	
		CSR A(a, 100, 100);

		for(size_t j = 0; j <= 100; j++){
			x0[j] = (std::rand()% (1000 - 0 + 1)) + 0;  // генерация начального приближения
		}
		
		for(size_t k = 0; k <= 100; k++){
			x[k] = (std::rand()% (1000 - 0 + 1)) + 0;  // генерация истинного решения
		}
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-15);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	std::vector<double> LU;
		
		std::size_t n = 0;
		auto begin = high_resolution_clock::now(); // начало записи времени
    		while (modul(res) > tolerance){
	    		LU = A.Jakobi_non_diagonal(y);
		    	y = mul(A.Jakobi_reverse_diag(), (b - LU));
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


void get_time_Gauss_Seidel(){
	std::srand(0);  //генерация случайных чисел
	double val; //значение из вектора a
	std::ofstream file;
        file.open("res_Gauss_Seidel.txt"); //файл для записи данных
	std::vector<double> a(10000); //вектор для генерации матрицы
	std::vector<double> x(100);
	std::vector<double> x0(100);
		for(size_t i = 0; i < 100; i++){
		   for(size_t j = 0; j < 100; j++){
			val = (std::rand()% (5000 - 10 + 1)) + 10; 
			if (i != j) a[i*100 + j] = val;
		   }
		}
		
		for(size_t i = 0; i < 100; i++){
		   double sum = 0.0;
		   for(size_t j = 0; j < 100; j++){
			if (i != j){
				sum += a[i*100 + j];
			}
		   }
		   a[i*100 + i] = sum + (std::rand()% (5000 - 10 + 1)) + 10;
		}
	
		CSR A(a, 100, 100);

		for(size_t j = 0; j <= 100; j++){
			x0[j] = (std::rand()% (1000 - 0 + 1)) + 0;  // генерация начального приближения
		}
		
		for(size_t k = 0; k <= 100; k++){
			x[k] = (std::rand()% (1000 - 0 + 1)) + 0;  // генерация истинного решения
		}
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-15);
		
		std::vector<double> y = x0; // x0 - начальное приближение
	    	std::vector<double> res = A * y - b;
	    	double diag; // диагональный элемент
		
		std::size_t n = 0;
		auto begin = high_resolution_clock::now(); // начало записи времени
    		while(modul(res) > tolerance){
		     for (std::size_t c = 0; c < y.size(); ++c) {
			 y[c] = b[c];
			 for(std::size_t k = A.get_rows()[c]; k < A.get_rows()[c+1]; k++){
			    if (A.get_cols()[k] == c){
			       diag = A(c,c);
			       continue;
			    }
			    
			    y[c] -= A.get_values()[k] * y[A.get_cols()[k]];
			 }
			 
			 y[c] /= diag;
		    }
		    
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



void get_time_Simple_iteration(){
	std::srand(0);  //генерация случайных чисел
	double val; //значение из вектора a
	std::ofstream file;
        file.open("res_Simple_iteration.txt"); //файл для записи данных
	std::vector<double> a(10000); //вектор для генерации матрицы
	std::vector<double> x(100);
	std::vector<double> x0(100);
		for(size_t i = 0; i < 100; i++){
		   for(size_t j = 0; j < 100; j++){
			val = (std::rand()% (5000 - 10 + 1)) + 10; 
			if (i != j) a[i*100 + j] = val;
		   }
		}
		
		for(size_t i = 0; i < 100; i++){
		   double sum = 0.0;
		   for(size_t j = 0; j < 100; j++){
			if (i != j){
				sum += a[i*100 + j];
			}
		   }
		   a[i*100 + i] = sum + (std::rand()% (5000 - 10 + 1)) + 10;
		}
	
		CSR A(a, 100, 100);

		for(size_t j = 0; j <= 100; j++){
			x0[j] = (std::rand()% (1000 - 0 + 1)) + 0;  // генерация начального приближения
		}
		
		for(size_t k = 0; k <= 100; k++){
			x[k] = (std::rand()% (1000 - 0 + 1)) + 0;  // генерация истинного решения
		}
	
		std::vector<double> b = A * x;  // вектор свободных членов
		double tolerance = pow(10,-4);
		double tau = 0.000000125;
		
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

int main(){
    get_time_Jacobi();
    get_time_Gauss_Seidel();
    get_time_Simple_iteration();
}
