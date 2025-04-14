#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>

TEST(Dense, TEST_1)
{
	
	std::vector<double> b = {1, 3, 2, 4, 7, 6, 2, 1, 0};
	std::vector<double> a = {1.1, 3.2, 2.3, 4.1, 7.2, 6.1, 2.0, 1.7, 0.8};
	Dense A(a, 3, 3);
	Dense B(b, 3, 3);
	Dense C = A + B;
	std::vector<double> d = {2.1, 6.2, 4.3, 8.1, 14.2, 12.1, 4.0, 2.7, 0.8};
	Dense D(d, 3, 3);
	for (std::size_t i = 0; i < 3; ++i) {
        	for (std::size_t j = 0; j < 3; ++j) {
       			ASSERT_NEAR(C(i,j), D(i,j) , 0.1);
        	}
    	}
}

TEST(Dense, TEST_2)
{
	std::vector<double> a = {0, 3, 2, 4, 7, 6, 2, 1, 0};
	Dense A(a, 3, 3);
	Dense B = A * 2.5;
	std::vector<double> d = {0, 7.5, 5.0, 10.0, 17.5, 15.0, 5.0, 2.5, 0.0};
	Dense D(d, 3, 3);
	for (std::size_t i = 0; i < 3; ++i) {
        	for (std::size_t j = 0; j < 3; ++j) {
       			ASSERT_NEAR(B(i,j), D(i,j) , 0.1);
        	}
    	}
}

TEST(Dense, TEST_3)
{
	std::vector<double> a = {1, 3, 2, 4, 7, 6, 2, 1, 0};
	std::vector<double> b = {1, 3, 2};
	Dense A(a, 3, 3);
	std::vector<double> c = A * b;
	std::vector<double> d = {14, 37, 5};
	
	for (std::size_t i = 0; i < 3; ++i) {
        	ASSERT_NEAR(c[i], d[i] , 0.1);
        }
}

TEST(Dense, TEST_4)
{
	std::vector<double> a = {2.5, 1, 3, 4.2};
	std::vector<double> b = {2, 1.4, 3.2, 1};
	Dense A(a, 2, 2);
	Dense B(b, 2, 2);
	Dense C = A * B;
	
	std::vector<double> d = {8.2, 4.5, 19.44, 8.4};
	Dense D(d, 2, 2);
	for (std::size_t i = 0; i < 2; ++i) {
        	for (std::size_t j = 0; j < 2; ++j) {
       			ASSERT_NEAR(C(i,j), D(i,j) , 0.1);
        	}
    	}
}
