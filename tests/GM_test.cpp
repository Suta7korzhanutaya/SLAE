#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>
#include <Gradient/GM.hpp>

TEST(Gradient_Descent, TEST_1)
{
	std::vector<double> a = {3, 1, 0, 2, 5, 1, 6, 0, 7};
	std::vector<double> x = {2.5, 1, 3};
	CSR A(a, 3, 3);
	double tolerance = pow(10,-13);
	std::vector<double> x0 = {-15, 100, 567};
	std::vector<double> b = A * x;
	
	std::vector<double> result = Gradient(A, b, x0, tolerance);	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}

TEST(Conjugate_Gradient, TEST_1)
{
	std::vector<double> a = {10, 16, 4, 17, 150, 30, 3, 30, 300};
	std::vector<double> x = {2.5, 12.5, 15.2};
	CSR A(a, 3, 3);
	
	double tolerance = pow(10,-13);
	std::vector<double> x0 = {13, 112, -74};
	std::vector<double> b = A * x;

	std::vector<double> result = Conjugate_Gradient(A, b, x0, tolerance);	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(result[i], x[i] , 0.1);
    	}
}
