#include "CSR.hpp"

CSR::CSR(const std::vector<double>& values, const std::vector<unsigned int>& cols, const std::vector<unsigned int>& rows) : values(values), cols(cols), rows(rows) {}


CSR::CSR(const std::vector<double>& data, std::size_t r, std::size_t c) {
	std::size_t temp = 0; // количество ненулевых элементов
	rows.resize(r+1);
	rows[0] = 0;
	for (std::size_t i = 0; i < r; i++) {
		for (std::size_t j = 0; j < c; j++) {
			if (data[i*c + j] != 0) {
				values.push_back(data[i*c + j]);
				cols.push_back(j);
				temp++;
			}
		}
		rows[i+1] = temp;
	}
}

double CSR::operator()(std::size_t i, std::size_t j) const {
    for (std::size_t k = rows[i]; k < rows[i + 1]; ++k) {
        if (cols[k] == j) {
            return values[k];
        }
    }
    return 0;
}

CSR CSR::operator*(double alpha) const {
    std::vector<double> result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = (values[i] * alpha);
    }
    CSR res(result, cols, rows);
    return res;
}

std::vector<double> CSR::operator*(const std::vector<double>& vec) const {
    std::vector<double> result(vec.size(), 0);
    for (std::size_t i = 0; i < vec.size(); ++i) {
        for (std::size_t k = rows[i]; k < rows[i + 1]; ++k) {
            result[i] += values[k] * vec[cols[k]];
        }
    }
    return result;
}

void CSR::print() const {
        std::cout << "values: ";
        for(std::size_t i = 0; i < values.size(); i++) std::cout << values[i] << " ";
        std::cout << std::endl;
        std::cout << "cols: ";
        for(std::size_t i = 0; i < cols.size(); i++) std::cout << cols[i] << " ";
        std::cout << std::endl;
        std::cout << "rows: ";
        for(std::size_t i = 0; i < rows.size(); i++) std::cout << rows[i] << " ";
        std::cout << std::endl;    
}

std::vector<double> CSR::Jakobi_reverse_diag() const {
	std::vector<double> v(rows.size()- 1);
	for (std::size_t i = 0; i < v.size(); ++i) {
        	for (std::size_t k = rows[i]; k < rows[i + 1]; ++k) {
        	    if (i == cols[k]){
            		v[i] = 1 / values[k];
            	    }
        	} 
	}
	return v;
}

std::vector<double> CSR::Jakobi_non_diagonal(const std::vector<double>& vec) const {
	std::vector<double> result(vec.size(), 0);
	for (std::size_t i = 0; i < vec.size(); ++i) {
        	for (std::size_t k = rows[i]; k < rows[i + 1]; ++k) {
        	    if (i == cols[k]) continue;
        	    else result[i] += values[k] * vec[cols[k]];        	    	
        	}
	}
	return result;
}

