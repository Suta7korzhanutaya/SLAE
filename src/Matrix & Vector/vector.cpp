#include "vector.hpp"

double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0;
    for (std::size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

std::vector<double> operator*(const std::vector<double>& a, const double& alpha) {
    std::vector<double> c(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
    	c[i] = alpha * a[i];
    }
    return c;
}

std::vector<double> operator*(const double& alpha, const std::vector<double>& a) {
    std::vector<double> c(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
    	c[i] = a[i] * alpha;
    }
    return c;
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> sum(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
    	sum[i] = a[i] + b[i];
    }
    return sum;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> res(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
    	res[i] = a[i] - b[i];
    }
    return res;
}

void print(const std::vector<double>& result) {
    for (std::size_t i = 0; i < result.size(); ++i) {
        std::cout << result[i] << " ";
        }
    std::cout << std::endl;    
}

std::vector<double> operator/(const std::vector<double>& a, const double& alpha) {
    std::vector<double> c(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
    	c[i] = a[i] / alpha ;
    }
    return c;
}



double modul(const std::vector<double>& x){
    double result = 0;
    for(std::size_t i = 0; i < x.size(); i++){
        result+= std::pow(x[i],2);
    }
    return std::pow(result,0.5);
}

double sq(const std::vector<double>& x){
    double result = 0;
    for(std::size_t i = 0; i < x.size(); i++){
        result+= std::pow(x[i],2);
    }
    return result;
}

int sign(const std::vector<double>& x){
    if(x[0] < 0) return -1;
    return 1;
}

std::vector<double> transposed(const std::vector<double>& x) {
    std::vector<double> result(x.size());
    for (std::size_t i = 0; i < x.size(); i++) {
        for (std::size_t j = 0; j < x.size(); j++) {
            result[i] = x[i] * x[j];
        }
    }
    return result;
}

std::vector<double> mul(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (std::size_t i = 0; i < a.size(); i++) {
            result[i] = a[i] * b[i];
        }
    return result;
}

