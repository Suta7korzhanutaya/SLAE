#include<gtest/gtest.h>
#include <Matrix & Vector/Dense.hpp>
#include <Solvers/Householder.cpp>
//#include <Solvers/1.cpp>
#include <cstdlib>

TEST(Householder, TEST_1)
{
	std::vector<double> a = {1, 0, 2, 3, 1, 7, 2, 4, 3};
	Dense A(a, 3, 3);
	Dense Q = Householder(A).first;
    	Dense R = Householder(A).second;
    	
	std::vector<double> q = {-0.267, -0.272, -0.925, 0.802, -0.469, 0.370, 0.535, 0.840, -0.092};
	std::vector<double> r = {-3.742, -2.940, 7.751, 0, 2.891, -1.310, 0, 0, 0.462};
	Dense B(q, 3, 3);
	Dense C(r, 3, 3);
	for (std::size_t i = 0; i < 3; ++i) {
        	for (std::size_t j = 0; j < 3; ++j) {
       			ASSERT_NEAR(Q(i,j), B(i,j) , 0.001);
       			ASSERT_NEAR(R(i,j), C(i,j) , 0.001);
        	}
    	}
}

TEST(Householder, TEST_2)
{
	std::vector<double> a = {1,2,4,0,0,5,0,3,6};
	Dense A(a, 3, 3);
	Dense Q = Householder(A).first;
    	Dense R = Householder(A).second;
    	Dense M = Q * R;
	
	EXPECT_TRUE(A == M);
}
