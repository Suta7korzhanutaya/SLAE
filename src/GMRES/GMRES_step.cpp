#include "GMRES.hpp"

std::pair<std::vector<double>, double> Step(const CSR& A, const std::vector<double>& b, std::vector<double> x0, const std::size_t m){

    Dense H(0, m + 1, m); //матрица Хессинберга
    Dense K(0, x0.size(), m + 1); //базис Крылова

    std::vector<std::pair<double, double>> rot(m);

    std::vector<double> r0 = A * x0 - b; //начальная невязка
    std::vector<double> e(m + 1, 0);
    e[0] = modul(r0);

    K.change_column(1/modul(r0) * r0, 0); //первый вектор Крылова - нормированная невязка

    for(std::size_t i = 0; i < m ; i++){
        Iteration(A, K, H, i);
        Givens_rot(H, rot, i);
       
        double e_m = rot[i].first * e[i] - rot[i].second * e[i+1];
        double e_p = rot[i].second * e[i] + rot[i].first * e[i+1];

        e[i] = e_m;
        e[i + 1] = e_p;   
    }

    std::vector<double> vec(e.begin(), e.begin() + m);
    
    std::vector<double> y = Back_Gauss(H, vec);
    std::vector<double> result(x0.size(), 0);
    
    for (std::size_t k = 0; k < y.size(); k++) {
    	result = result + y[k] * K.get_column(k); //ЛК базиса Крылова
    }
    
    double ans = (e[m] < 0) ? -e[m] : e[m]; //абсолютное значение элемента
    return {x0 - result, ans};
}
