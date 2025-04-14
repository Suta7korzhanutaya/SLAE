#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>
#include <SIM/SIM.hpp>

TEST(Simple_iteration, TEST_1)
{
	std::vector<double> a = {3, 1, 0, 2, 5, 1, 6, 0, 7};
	std::vector<double> x = {2.5, 1, 3};
	CSR A(a, 3, 3);
	double tolerance = pow(10,-13);
	double tau = 0.111;
	std::vector<double> x0 = {-15, 100, 567};
	std::vector<double> b = A * x;
	
	std::vector<double> result = Simple(A, b, x0, tau, tolerance);	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Gauss_Seidel, TEST_1)
{
	std::vector<double> a = {33.7, 2, 11, 20, 42, 3, 1, 0, 7};
	std::vector<double> x = {0.47, 5.24, 3.18};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	std::vector<double> x0 = {-24, 5.678, 517};
	
	std::vector<double> b = A * x;
	std::vector<double> result = Gauss_Seidel(A, b, x0, tolerance);
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Jacobi, TEST_1)
{
	std::vector<double> a = {4, 1, 1, 0, 5, 2, 1, 4, 7};
	std::vector<double> x = {1, 1, 1};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	std::vector<double> x0 = {2, 2, 2};
	
	std::vector<double> b = A * x;
	std::vector<double> result = Jacobi(A, b, x0, tolerance);
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Chebyshev, TEST_1)
{
	std::vector<double> a = {3, 1, 0, 2, 5, 1, 6, 0, 7};
	std::vector<double> x = {2.5, 1, 3};
	CSR A(a, 3, 3);
	double tolerance = pow(10,-12);

	std::vector<double> x0 = {-15, 260, 27};
	std::vector<double> b = A * x;
	
	std::vector<double> result = Chebyshev(A, b, x0, tolerance, 8, 2.738, 7.602);
		
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Sym_Gauss_Seidel, TEST_1)
{
	std::vector<double> a = {33.7, 2, 11, 20, 42, 3, 1, 0, 7};
	std::vector<double> x = {0.47, 5.24, 3.18};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	std::vector<double> x0 = {-24, 5.678, 517};
	
	std::vector<double> b = A * x;
	std::vector<double> result = Sym_Gauss_Seidel(A, b, x0, tolerance);
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Chebyshev_Sym_GS, TEST_1)
{
	std::vector<double> a = {12.5, 1, 3, 9.7, 52, 11.2, 0, 8.7, 67};
	std::vector<double> x = {0.4, 5.2, 3.7};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	double rho = 0.08;
	std::vector<double> x0 = {-24, 5, 17};
	
	std::vector<double> b = A * x;
	std::vector<double> result = Chebyshev_Sym_GS(A, b, x0, tolerance, rho);
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Chebyshev_Sym_GS, TEST_2)
{
	std::vector<double> a = {21, 5, 3, 3, 16, 3, 5, 1, 19};
	std::vector<double> x = {0.009, 0.288, 0.456};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-12);
	double rho = 0.073381;
	std::vector<double> x0 = {0, 0, 0};
	
	std::vector<double> b = {3, 6, 9};
	std::vector<double> result = Chebyshev_Sym_GS(A, b, x0, tolerance, rho);
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

