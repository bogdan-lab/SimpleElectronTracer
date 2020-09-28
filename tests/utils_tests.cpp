#include <gtest/gtest.h>
#include <cmath>
#include "utils.hpp"


TEST(UtilsTests, ConstructorTests){
    Vector v;
    EXPECT_EQ(v.GetX(), 0.0);
    EXPECT_EQ(v.GetY(), 0.0);
    EXPECT_EQ(v.GetZ(), 0.0);
    Vector v1(1.0, 2.0, 3.0);
    EXPECT_EQ(v1.GetX(), 1.0);
    EXPECT_EQ(v1.GetY(), 2.0);
    EXPECT_EQ(v1.GetZ(), 3.0);
    Vector v2(Vector(4.0, 8.0, 12.0), Vector(1.0, 3.0, 5.0));
    EXPECT_EQ(v2.GetX(), -3.0);
    EXPECT_EQ(v2.GetY(), -5.0);
    EXPECT_EQ(v2.GetZ(), -7.0);
}


TEST(UtilsTests, DotTest){
    Vector v1(2.0, 5.0, 9.0);
    Vector v2(1.0, 3.0, -0.5);
    EXPECT_EQ(v1.Dot(v2), 12.5);
    EXPECT_EQ(v2.Dot(v1), 12.5);
}

TEST(UtilsTests, CrossTest){
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(4.0, 5.0, 6.0);
    Vector v3 = v1.Cross(v2);
    EXPECT_EQ(v3.GetX(), -3.0);
    EXPECT_EQ(v3.GetY(), 6.0);
    EXPECT_EQ(v3.GetZ(), -3.0);
    Vector v4 = v2.Cross(v1);
    EXPECT_EQ(v4.GetX(), 3.0);
    EXPECT_EQ(v4.GetY(), -6.0);
    EXPECT_EQ(v4.GetZ(), 3.0);
}

TEST(UtilsTests, TimesTest){
    Vector v(1.0, 2.0, 3.0);
    Vector res = v.Times(-5);
    EXPECT_EQ(res.GetX(), -5.0);
    EXPECT_EQ(res.GetY(), -10.0);
    EXPECT_EQ(res.GetZ(), -15.0);
}


TEST(UtilsTests, NormTest){
    Vector v(1.0, 2.0, -3.0);
    Vector res = v.Norm();
    EXPECT_EQ(res.GetX(), 1.0/std::sqrt(14.0));
    EXPECT_EQ(res.GetY(), 2.0/std::sqrt(14.0));
    EXPECT_EQ(res.GetZ(), -3.0/std::sqrt(14.0));
}

int main(int argc, char* argv[]){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

