#include <iostream>
#include <vector>

std::vector<double> solveTridiagonalMatrix(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d)
{
    int n = d.size();
    std::vector<double> p(n);
    std::vector<double> q(n);
    std::vector<double> x(n);

    int N = n-1; 

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for (int i = 1; i < N; ++i)
    {
        p[i] = - c[i] / (a[i]*p[i-1] + b[i]);
        q[i] = (d[i] - a[i] * q[i-1]) / (a[i]*p[i-1] + b[i]);
    }

    x[N] = (d[N] - a[N] * q[N-1]) / (a[N]*p[N-1] + b[N]);

    for (int i = N-1 ; i >= 0; --i)
    {
        x[i] = p[i]*x[i+1] + q[i];
    }

    return x;
}

int main()
{
    int n;
    std::cout << "Size of TridiagonalMatrix: ";
    std::cin >> n;

    std::vector<double> a(n); 
    std::cout << "Lower diagonal:" << std::endl; 
    for (int i = 1; i < n; ++i) {
        std::cin >> a[i];
    }
    
    std::vector<double> b(n);
    std::cout << "Main diagonal:" << std::endl; 
    for (int i = 0; i < n; ++i) {
        std::cin >> b[i];
    }

    std::vector<double> c(n); 
    std::cout << "Upper diagonal:" << std::endl; 
    for (int i = 0; i < n-1; ++i) {
        std::cin >> c[i];
    }

    std::vector<double> d(n); 
    std::cout << "Free members:" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cin >> d[i];
    }

    std::vector<double> solution = solveTridiagonalMatrix(a, b, c, d);

    std::cout << "Solution:" << std::endl; 
    for (int i = 0; i < solution.size(); ++i) {
        std::cout << "x[" << i << "] = " << solution[i] << std::endl;
    }

    return 0;
}
