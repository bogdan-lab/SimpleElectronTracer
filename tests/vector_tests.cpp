#include <gtest/gtest.h>
#include <cmath>
#include "utils.hpp"


TEST(Vec3Tests, ConstructorTests){
    Vec3 v;
    EXPECT_EQ(v.GetX(), 0.0);
    EXPECT_EQ(v.GetY(), 0.0);
    EXPECT_EQ(v.GetZ(), 0.0);
    Vec3 v1(1.0, 2.0, 3.0);
    EXPECT_EQ(v1.GetX(), 1.0);
    EXPECT_EQ(v1.GetY(), 2.0);
    EXPECT_EQ(v1.GetZ(), 3.0);
    Vec3 v2(Vec3(4.0, 8.0, 12.0), Vec3(1.0, 3.0, 5.0));
    EXPECT_EQ(v2.GetX(), -3.0);
    EXPECT_EQ(v2.GetY(), -5.0);
    EXPECT_EQ(v2.GetZ(), -7.0);
}


TEST(Vec3Tests, DotTest){
    Vec3 v1(2.0, 5.0, 9.0);
    Vec3 v2(1.0, 3.0, -0.5);
    EXPECT_EQ(v1.Dot(v2), 12.5);
    EXPECT_EQ(v2.Dot(v1), 12.5);
}

TEST(Vec3Tests, CrossTest){
    Vec3 v1(1.0, 2.0, 3.0);
    Vec3 v2(4.0, 5.0, 6.0);
    Vec3 v3 = v1.Cross(v2);
    EXPECT_EQ(v3.GetX(), -3.0);
    EXPECT_EQ(v3.GetY(), 6.0);
    EXPECT_EQ(v3.GetZ(), -3.0);
    Vec3 v4 = v2.Cross(v1);
    EXPECT_EQ(v4.GetX(), 3.0);
    EXPECT_EQ(v4.GetY(), -6.0);
    EXPECT_EQ(v4.GetZ(), 3.0);
}

TEST(Vec3Tests, TimesTest){
    Vec3 v(1.0, 2.0, 3.0);
    Vec3 res = v.Times(-5);
    EXPECT_EQ(res.GetX(), -5.0);
    EXPECT_EQ(res.GetY(), -10.0);
    EXPECT_EQ(res.GetZ(), -15.0);
}


TEST(Vec3Tests, GetDistanceTest){
    Vec3 v(1.0, 2.0, 3.0);
    Vec3 v2(1.0, 6.0, 0.0);
    EXPECT_EQ(v2.GetDistance(v), 5.0);
    EXPECT_EQ(v2.GetDistance(v), v.GetDistance(v2));
}

TEST(Vec3Tests, NormTest){
    Vec3 res(1.0, 2.0, -3.0);
    res.Norm();
    EXPECT_EQ(res.GetX(), 1.0/std::sqrt(14.0));
    EXPECT_EQ(res.GetY(), 2.0/std::sqrt(14.0));
    EXPECT_EQ(res.GetZ(), -3.0/std::sqrt(14.0));
}

TEST(Vec3Tests, TestCompare){
    Vec3 lhs(2.0, 15.0, -9);
    Vec3 rhs1(2.0, 15.0, -9);
    Vec3 rhs2(-2.0, 15.0, 9);
    EXPECT_TRUE(lhs==rhs1);
    EXPECT_FALSE(lhs==rhs2);
}


TEST(Vec3Tests, Arithmetic){
    Vec3 lhs(1.0, 2.0, 3.0);
    Vec3 rhs(25.0, 17.0, 10.0);
    EXPECT_TRUE(lhs+rhs == Vec3(26.0, 19.0, 13.0));
    EXPECT_TRUE(lhs+rhs == rhs+lhs);
    EXPECT_TRUE(rhs-lhs == Vec3(24.0, 15.0, 7.0));
    EXPECT_TRUE(lhs-rhs == (rhs-lhs).Times(-1));
}
