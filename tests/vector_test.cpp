#include<gtest/gtest.h>
#include <Matrix & Vector/vector.hpp>

TEST(vector, TEST_1)
{
	std::vector<double> a = {11.1, 18.7, 5.68};
	std::vector<double> b = {9, 1, 5};
	std::vector<double> c = {7, 3, 13};
	double d = (2 * a) * (c + b);
	ASSERT_NEAR(d, 709.28 , 0.1);
}

TEST(vector, TEST_2)
{
	std::vector<double> a = {3, 4, 5};
	std::vector<double> b = {2.2, 3.4, 5.1};
	std::vector<double> c = {7, 3, 4};
	double d = (b * 2) * (c - a);
	ASSERT_NEAR(d, 0.6 , 0.1);
}

TEST(vector, TEST_3)
{
	std::vector<double> a = {7.2, 3.4, 8.1, 7.0};
	std::vector<double> b = {2.7, 9.1, 2.25, 3.11};
	std::vector<double> d = (3 * b) + (a * 2);
	ASSERT_NEAR(d[0], 22.5, 0.1);
    	ASSERT_NEAR(d[1], 34.1, 0.1);
    	ASSERT_NEAR(d[2], 22.95, 0.01);
    	ASSERT_NEAR(d[3], 23.33, 0.01);
}
