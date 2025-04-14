#include "Dense.hpp"

Dense::Dense(const double elem, const std::size_t rows, const std::size_t cols):
          data(std::vector(rows * cols, elem)), rows(rows), cols(cols) {}

Dense::Dense(const std::vector<double>& data, std::size_t rows, std::size_t cols) : rows(rows), cols(cols), data(data) {}

double Dense::operator()(std::size_t i, std::size_t j) const {
    return data[i*cols + j];
}

std::size_t Dense::get_cols() const {
	return cols;
}
	    
std::size_t Dense::get_rows() const {
	return rows;
}


Dense Dense::operator+(const Dense& other) const {
    std::vector<double> result(data.size(), 0); 
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result[i* cols + j] = data[i * cols + j] + other(i,j);
        }
    }
    Dense ans(result, cols, rows);
    return ans;
}

Dense Dense::operator*(double alpha) const {
    std::vector<double> result(data.size(), 0);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result[i* cols + j] = data[i * cols + j] * alpha;
        }
    }
    Dense ans(result, cols, rows);
    return ans;
}

std::vector<double> Dense::operator*(const std::vector<double>& vec) const {
    std::vector<double> result(rows);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            result[i] += data[i * cols + j] * vec[j];
        }
    } 
    return result;
}

Dense Dense::operator*(const Dense& other) const {
    std::vector<double> result(data.size(), 0);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < other.cols; ++j) {
            for (std::size_t k = 0; k < cols; ++k) {
                result[i * cols + j] += data[i * cols + k] * other(k, j);
            }
        }
    }
    Dense ans(result, cols, rows);
    return ans;
}

void Dense::print() const {
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            std::cout << data[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<double> Dense::get_elements(const std::vector<double>& data) const {
    return data;
}

const std::vector<double> Dense::get_column(std::size_t num) const{
        std::vector<double> column(rows, 0);
        for(int i = 0; i < rows; i++){
            column[i] = data[cols * i + num];
        }
        return column;
    }

const std::vector<double> Dense::get_row(std::size_t num) const{
        std::vector<double> row(cols, 0);
        for(int j = 0; j < cols; j++){
            row[j] = data[num*cols + j];
        }
        return row;
}


double& Dense::get_element(std::size_t i,  std::size_t j){
        return(data[i * cols + j]);
    }
    
   

void Dense::change_column(const std::vector<double>& vec, std::size_t i){
	for(std::size_t j = 0; j < rows; j++){
		data[cols* j + i] = vec[j];
        }
}



/*void Dense::single(std::size_t n) {
    std::vector<double> vec(n * n);
    for (std::size_t i = 0; i < n; i++) {
        vec[n * i + i] = 1;
    }
    return Dense A(vec, n, n);
}
*/

void Dense::swap_elements(std::size_t i, std::size_t j, double num) {
        data[i * cols + j] = num;
    }

Dense Dense::transpon() const{
    std::vector<double> v(cols*rows);
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 1 + i; j < cols; j++) {
            v[j * cols + i] = data[i * cols + j];
            v[i * cols + j] = data[j * cols + i];
        }
        v[i * cols + i] = data[i * cols + i];
    }
    Dense N(v, rows, cols);
    return N;
}

bool Dense::operator == (const Dense & other){
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
        if (std::abs(data[i * cols + i] - other(i, j)) > 0.01)
            return false;
        }
    }
    return true;
}
