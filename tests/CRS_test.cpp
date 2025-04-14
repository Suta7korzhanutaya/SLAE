#include<gtest/gtest.h>
#include <Matrix & Vector/CSR.hpp>

TEST(CSR, TEST_1)
{
	std::vector<double> a = {1, 0, 2, 0, 5, 4, 2, 1, 0};
	std::vector<double> v = {1, 2, 3, 4, 5, 1};
	std::vector<unsigned int> c = {0, 0, 1, 2, 0, 2};
	std::vector<unsigned int> r = {0, 1, 4, 6};
	CSR D(a, 3, 3);
	CSR B(v, c, r);
	
	CSR T = D * 2;
	
	std::vector<double> t = {2, 1, 2};
	std::vector<double> b = T * t;
	
	std::vector<double> ans = {12, 26, 10};
	for (std::size_t i = 0; i < 3; i ++) {
        	ASSERT_NEAR(b[i], ans[i] , 0);
    	}
    	
    	std::vector<double> vec = {1 , 2, 3};
	std::vector<double> d = B * vec;
	std::vector<double> res = {1, 20, 8};
	for (std::size_t i = 0; i < res.size(); i ++) {
        	ASSERT_NEAR(d[i], res[i] , 0.1);
    	}
}
