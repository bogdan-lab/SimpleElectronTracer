#include <gtest/gtest.h>
#include "math.hpp"


TEST(MatrixTests, GenerationTest){
    Vec3 i(0.1, 0.0, 0.0);
    Vec3 j(0.0, 25.0, 0.0);
    Vec3 k(0.0, 0.0, -15.0);
    ONBasis_3x3 m(i, j, k);
    EXPECT_TRUE(m.GetXVec()==i.Norm());
    EXPECT_TRUE(m.GetYVec()==j.Norm());
    EXPECT_TRUE(m.GetZVec()==k.Norm());
    auto gen_test = [](const Vec3& i, const Vec3& j, const Vec3& k){
                            return ONBasis_3x3(i, j, k);
                        };
    Vec3 k_wrong(1.0, 2.0, 3.0);
    EXPECT_DEATH(gen_test(i, j, k_wrong), "Basis is not orthogonal!");

    Vec3 x(1.0, 2.0, 3.0);
    Vec3 y=x.Cross(Vec3(1.0, 0.0, 0.0));
    Vec3 z = x.Cross(y);
    ONBasis_3x3 m2(x,y,z);
    EXPECT_TRUE(m2.GetXVec()==x.Norm());
    EXPECT_TRUE(m2.GetYVec()==y.Norm());
    EXPECT_TRUE(m2.GetZVec()==z.Norm());
}


TEST(MatrixTests, CoordinateTransitionTest){
    ONBasis_3x3 m({1.0, 0.0, 0.0},
                  {0.0, 1.0, 0.0},
                  {0.0, 0.0, -1.0});
    Vec3 tst(0.1, 0.2, 0.3);
    Vec3 res = m.ApplyToVec(tst);
    EXPECT_EQ(res.GetX(), tst.GetX());
    EXPECT_EQ(res.GetY(), tst.GetY());
    EXPECT_EQ(res.GetZ(), -tst.GetZ());
}


TEST(MatrixTests, TestGenerationFromZ){
    Vec3 new_z(0.1, 0.2, 0.3);
    ONBasis_3x3 m(new_z);
    Vec3 x = m.GetXVec();
    Vec3 y = m.GetYVec();
    Vec3 z = m.GetZVec();

    EXPECT_NEAR(0.0, x.Dot(y), 1e-15);
    EXPECT_NEAR(0.0, y.Dot(z), 1e-15);
    EXPECT_NEAR(0.0, x.Dot(z), 1e-15);

    Vec3 z_orig(0.0, 0.0, 1.0);
    Vec3 z_trans = m.ApplyToVec(z_orig);
    new_z = new_z.Norm();
    //Check if normed basis generated according to Zvec will transform
    //original z vector into itself
    EXPECT_NEAR(new_z.GetX(), z_trans.GetX(), 1e-15);
    EXPECT_NEAR(new_z.GetY(), z_trans.GetY(), 1e-15);
    EXPECT_NEAR(new_z.GetZ(), z_trans.GetZ(), 1e-15);
}

TEST(MatrixTests, FromOriginalCoorsToThisAndBack){
    Vec3 new_x(-1.0, 0.0, 0.0);
    Vec3 new_y(0.0, -1.0, 0.0);
    Vec3 new_z(0.0, 0.0, -1.0);
    ONBasis_3x3 m(new_x, new_y, new_z);
    Vec3 tst(1.0, 2.0, 3.0);
    Vec3 tst_basis = m.FromOriginalCoorsToThis(tst);
    EXPECT_EQ(tst_basis.GetX(), -tst.GetX());
    EXPECT_EQ(tst_basis.GetY(), -tst.GetY());
    EXPECT_EQ(tst_basis.GetZ(), -tst.GetZ());
    Vec3 tst_back = m.FromThisCoorsToOriginal(tst_basis);
    EXPECT_EQ(tst_back.GetX(), tst.GetX());
    EXPECT_EQ(tst_back.GetY(), tst.GetY());
    EXPECT_EQ(tst_back.GetZ(), tst.GetZ());
}
