#include "GMRES.hpp"

void Iteration(const CSR& A, Dense& K, Dense& H, const std::size_t i){
    std::vector<double> v = A * K.get_column(i);
    double h;
    
    //заполняем матрицу Хессинберга H
    for (std::size_t k = 0; k <= i; k++){
        h = v * K.get_column(k);
        v = v - h * K.get_column(k);
        H.get_element(k, i) = h;
    }

    h = modul(v);
    H.get_element(i + 1, i) = h;
    v = 1 / h * v;
    K.change_column(v, i + 1);
}


void Givens_rot(Dense& H, std::vector<std::pair<double, double>>& rot, const std::size_t i){
    //применяем все предыдущие вращения к новому столбцу H
    for (std::size_t k = 0; k < i; k++){
        double h = rot[k].first * H(k, i) - rot[k].second * H(k + 1, i);
        double h_n = rot[k].first * H(k+1, i) + rot[k].second * H(k, i);
        H.get_element(k, i) = h;
        H.get_element(k + 1, i) = h_n;
    }
    
    
    //новое вращение к столбцу H
    double norm = std::pow(H(i, i)* H(i, i) + H(i + 1, i) * H(i+1, i) , 0.5);
    double cos_th = H(i, i) / norm;
    double sin_th = -H(i + 1, i) / norm;
  
    double h  = cos_th * H(i, i) - sin_th * H(i + 1, i);
    double h_n = sin_th * H(i, i) + cos_th * H(i + 1, i);

    H.get_element(i, i) = h;
    H.get_element(i + 1, i) = h_n;

    rot[i] = (std::make_pair(cos_th, sin_th));
}

