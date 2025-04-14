#include <gtest/gtest.h>
#include "solve.hpp"


TEST(TRIDIAGONAL, TEST_1) 
{
    std::vector<double> a{6.7, 0.75};
    std::vector<double> b{10.5, 9.4, 3.4};
    std::vector<double> c{7.8, 5.2};
    std::vector<double> d{11.1, 18.7, 5.68};
    
    std::vector<double> x = solveTridiagonalMatrix(a, b, c, d);
    ASSERT_NEAR(x[0], 0.3926838028, 1e-10);
    ASSERT_NEAR(x[1], 0.8944641116, 1e-10);
    ASSERT_NEAR(x[2], 1.4732799754, 1e-10);
}

TEST(TRIDIAGONAL, TEST_2) 
{
    std::vector<double> a{9, 1, 5};
    std::vector<double> b{42, 3, 2, 24};
    std::vector<double> c{87, 4, 6};
    std::vector<double> d{7, 8, 6, 8};
    
    std::vector<double> x = solveTridiagonalMatrix(a, b, c, d);
    ASSERT_NEAR(x[0], -1.298146046, 1e-7);
    ASSERT_NEAR(x[1], 0.7071509648, 1e-7);
    ASSERT_NEAR(x[2], 4.39046538, 1e-7);
    ASSERT_NEAR(x[3], -0.5813469542, 1e-7);
}

TEST(TRIDIAGONAL, TEST_3) 
{
    std::vector<double> a{18, 9, 28, 19};
    std::vector<double> b{24, 57, 12, 68, 37};
    std::vector<double> c{7, 3, 13, 43};
    std::vector<double> d{13, 78, 63, 14, 9};
    
    std::vector<double> x = solveTridiagonalMatrix(a, b, c, d);
    ASSERT_NEAR(x[0], 0.3956750001, 1e-10);
    ASSERT_NEAR(x[1], 0.5005428569, 1e-10);
    ASSERT_NEAR(x[2], 14.11563572, 1e-7);
    ASSERT_NEAR(x[3], -8.530193411, 1e-9);
    ASSERT_NEAR(x[4], 4.623612833, 1e-9);
}
