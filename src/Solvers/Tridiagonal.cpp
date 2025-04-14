#include "Tridiagonal.hpp"

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
			p[i] = - c[i] / (a[i-1]*p[i-1] + b[i]);
			q[i] = (d[i] - a[i-1] * q[i-1]) / (a[i-1]*p[i-1] + b[i]);
	    	}

	    	x[N] = (d[N] - a[N-1] * q[N-1]) / (a[N-1]*p[N-1] + b[N]);
		
	    	for (int i = N-1 ; i >= 0; --i)
	    	{
			x[i] = p[i]*x[i+1] + q[i];
	    	}
		return x;
};
