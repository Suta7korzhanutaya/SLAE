#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>
#include <Solvers/Gauss.cpp>

TEST(Gauss, TEST_1)
{
	std::vector<double> b = {19, 19, 14};
	std::vector<double> a = {1, 3, 4, 2, 7, 1, 5, 0, 3};
	Dense A(a, 3, 3);
	
	std::vector<double> x = Gauss(A, b);
	std::vector<double> res = {1, 2, 3};
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(x[i], res[i] , 0.1);
    	}
}

TEST(Gauss, TEST_2)
{
	std::vector<double> b = {30, 8, 42};
	std::vector<double> a = {0.5, 1, 2, 1, 0.5, 0, 2, 3, 1};
	Dense A(a, 3, 3);
	
	std::vector<double> x = Gauss(A, b);
	std::vector<double> res = {4, 8, 10};
	
	for (std::size_t i = 0; i < 3; ++i) {
       		ASSERT_NEAR(x[i], res[i] , 0.1);
    	}
}
