#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include "utils.hpp"
#include "particle.hpp"
#include "reflector.hpp"

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
    Vec3 v(1.0, 2.0, -3.0);
    Vec3 res = v.Norm();
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

TEST(MatrixTests, GenerationTest){
    Vec3 i(0.2, 0.5, 20.0);
    Vec3 j(0.5, 6.5, 20.0);
    Vec3 k(0.24, 0.5, 2.0);
    Basis_3x3 m(i, j, k);
    std::vector<Vec3> basis = m.GetBasisCols();
    EXPECT_TRUE(basis[0]==i);
    EXPECT_TRUE(basis[1]==j);
    EXPECT_TRUE(basis[2]==k);
}


TEST(MatrixTests, NormalizationTest){
    Vec3 i(0.2, 0.5, 20.0);
    Vec3 j(0.5, 6.5, 20.0);
    Vec3 k(0.24, 0.5, 2.0);
    Basis_3x3 m(i, j, k);
    m = m.Norm();
    std::vector<Vec3> basis = m.GetBasisCols();
    EXPECT_TRUE(basis[0]==i.Norm());
    EXPECT_TRUE(basis[1]==j.Norm());
    EXPECT_TRUE(basis[2]==k.Norm());
}

TEST(MatrixTests, DeterminantTest){
     Basis_3x3 m({0.6, 0.8, 0.0},
                  {0.0, 0.6, -0.8},
                  {-1.0, 0.0, 0.0});
    EXPECT_NEAR(0.64, m.GetDeterminant(), 1e-15);
}

TEST(MatrixTests, InverseTest){
    Basis_3x3 m({0.2, 0.8, 0.0},
                {0.0, 4.0, 4.0},
                {-0.5, 0.0, 1.0});
    Basis_3x3 inv = m.GetInverse();
    std::vector<Vec3> inv_basis = inv.GetBasisCols();
    EXPECT_LE(inv_basis[0].GetDistance(Vec3(-5, 1, -4)), 1e-15)
            << inv_basis[0];
    EXPECT_LE(inv_basis[1].GetDistance(Vec3(2.5, -0.25, 1)), 1e-15)
            << inv_basis[1];
    EXPECT_LE(inv_basis[2].GetDistance(Vec3(-2.5, 0.5, -1)), 1e-15)
            << inv_basis[2];
    Basis_3x3 m2({2.0, -2.0, 1.0},
                  {2.0, 1.0, -2.0},
                  {2.0, 1.0, -2.0});
    ASSERT_DEATH(m2.GetInverse(), "Basis have zero deternminant!");
}

TEST(MatrixTests, TransposeTest){
    Basis_3x3 m({0.6, 0.8, 0.0},
                  {0.0, 0.6, -0.8},
                  {-1.0, 0.0, 0.0});
    Basis_3x3 tr = m.Transpose();
    std::vector<Vec3> tr_basis = tr.GetBasisCols();
    EXPECT_TRUE(tr_basis[0]==Vec3(0.6, 0.0, -1.0)) << tr_basis[0];
    EXPECT_TRUE(tr_basis[1]==Vec3(0.8, 0.6, 0.0)) << tr_basis[1];
    EXPECT_TRUE(tr_basis[2]==Vec3(0.0, -0.8, 0.0)) << tr_basis[2];
}


TEST(MatrixTests, CoordinateTransitionTest){
    Basis_3x3 m({1.0, 0.0, 0.0},
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
    Basis_3x3 m(new_z);
    std::vector<Vec3> basis = m.GetBasisCols();

    EXPECT_NEAR(0.0, basis[0].Dot(basis[1]), 1e-15);
    EXPECT_NEAR(0.0, basis[1].Dot(basis[2]), 1e-15);
    EXPECT_NEAR(0.0, basis[0].Dot(basis[2]), 1e-15);

    Vec3 z_orig(0.0, 0.0, 1.0);
    m = m.Norm();
    Vec3 z_trans = m.ApplyToVec(z_orig);
    new_z = new_z.Norm();
    //Check if normed basis generated according to Zvec will transform
    //original z vector into itself
    EXPECT_NEAR(new_z.GetX(), z_trans.GetX(), 1e-15);
    EXPECT_NEAR(new_z.GetY(), z_trans.GetY(), 1e-15);
    EXPECT_NEAR(new_z.GetZ(), z_trans.GetZ(), 1e-15);
}


TEST(ParticleTests, ParticleGenerationTest1){
        Particle pt;
        EXPECT_EQ(pt.GetPosition().GetX(), 0.0);
        EXPECT_EQ(pt.GetPosition().GetY(), 0.0);
        EXPECT_EQ(pt.GetPosition().GetZ(), 0.0);
        EXPECT_EQ(pt.GetDirection().GetX(), 0.0);
        EXPECT_EQ(pt.GetDirection().GetY(), 0.0);
        EXPECT_EQ(pt.GetDirection().GetZ(), 0.0);
        EXPECT_EQ(pt.GetSurfCount(), 0);
        EXPECT_EQ(pt.GetVolCount(), 0);
    }

TEST(ParticleTests, ParticleGenerationTest2){
        Particle pt(Vec3(1.0, 2.0, 3.0), Vec3(4.0, 5.0, 6.0));
        EXPECT_EQ(pt.GetPosition().GetX(), 1.0);
        EXPECT_EQ(pt.GetPosition().GetY(), 2.0);
        EXPECT_EQ(pt.GetPosition().GetZ(), 3.0);
        EXPECT_NEAR(pt.GetDirection().GetX(), 4.0/sqrt(77.0), 1e-15);
        EXPECT_NEAR(pt.GetDirection().GetY(), 5.0/sqrt(77.0), 1e-15);
        EXPECT_NEAR(pt.GetDirection().GetZ(), 6.0/sqrt(77.0), 1e-15);
        EXPECT_EQ(pt.GetSurfCount(), 0);
        EXPECT_EQ(pt.GetVolCount(), 0);
    }

TEST(ParticleTests, ParticleGenerationTest3){
        std::mt19937 rnd_gen(42);
        Vec3 pos(-1.0, 15.0, 48.0);
        Vec3 dir(-2.0, 5.0, 10.0);
        Particle ptst(pos, dir, rnd_gen);
        EXPECT_EQ(ptst.GetPosition().GetX(), -1.0);
        EXPECT_EQ(ptst.GetPosition().GetY(), 15.0);
        EXPECT_EQ(ptst.GetPosition().GetZ(), 48.0);
        EXPECT_EQ(ptst.GetSurfCount(), 0);
        EXPECT_EQ(ptst.GetVolCount(), 0);
        for(size_t i=0; i<100; i++){
            Particle pt(pos, dir, rnd_gen);
            EXPECT_GE(dir.Dot(pt.GetDirection()), 0);
            EXPECT_NEAR(pt.GetDirection().Length(), 1.0, 1e-15);
        }
}

TEST(ParticleTests, GetRandomVelTest){
    std::mt19937 rnd_gen(42);
    Vec3 dir(1.0, -5.0, 8.0);
    Particle pt;
    for(size_t i=0; i<100; i++){
        Vec3 res = pt.GetRandomVel(dir, rnd_gen);
        EXPECT_GE(dir.Dot(res), 0);
        EXPECT_NEAR(res.Length(), 1.0, 1e-15);
    }
}

TEST(ParticleTests, MakeGasCollisionTest){
    double distance = 15;
    Vec3 dir(4.0, 5.0, 6.0);
    Vec3 pos(1.0, 2.0, 3.0);
    Particle pt(pos, dir);
    std::mt19937 rnd_gen(42);
    pt.MakeGasCollision(distance, rnd_gen);
    EXPECT_EQ(pt.GetVolCount(), 1);
    EXPECT_EQ(pt.GetSurfCount(), 0);
    dir = dir.Norm();
    EXPECT_NEAR(pt.GetPosition().GetX(), pos.GetX() + distance*dir.GetX(), 1e-15);
    EXPECT_NEAR(pt.GetPosition().GetY(), pos.GetY() + distance*dir.GetY(), 1e-15);
    EXPECT_NEAR(pt.GetPosition().GetZ(), pos.GetZ() + distance*dir.GetZ(), 1e-15);
    EXPECT_NE(dir.GetX(), pt.GetDirection().GetX());
    EXPECT_NE(dir.GetY(), pt.GetDirection().GetY());
    EXPECT_NE(dir.GetZ(), pt.GetDirection().GetZ());
}


TEST(ReflectorTests, MirrorReflectorTest){
    Vec3 normal(0.0, 0.0, -1.0);
    Vec3 dir(2.0, 3.0, 5.0);
    Particle pt({0.1, 0.2, 0.3}, dir);
    std::mt19937 rnd_gen(42);
    {
        MirrorReflector test(0.0);
        auto res = test.ReflectParticle(pt, normal, rnd_gen);
        EXPECT_FALSE(res.has_value());
    }
    {
        MirrorReflector test(1.0);
        auto res = test.ReflectParticle(pt, normal, rnd_gen);
        EXPECT_TRUE(res.has_value());
        dir = dir.Norm();
        EXPECT_NEAR(res->GetX(), dir.GetX(), 1e-15);
        EXPECT_NEAR(res->GetY(), dir.GetY(), 1e-15);
        EXPECT_NEAR(res->GetZ(), -dir.GetZ(), 1e-15);
    }
}


TEST(ReflectorTests, LambertianReflectorTest){
    Vec3 normal(0.0, 0.0, -1.0);
    Vec3 dir(2.0, 3.0, 5.0);
    Particle pt({0.1, 0.2, 0.3}, dir);
    std::mt19937 rnd_gen(42);
    {
        LambertianReflector test(0.0);
        auto res = test.ReflectParticle(pt, normal, rnd_gen);
        EXPECT_FALSE(res.has_value());
    }
    {
        LambertianReflector test(1.0);
        auto res = test.ReflectParticle(pt, normal, rnd_gen);
        EXPECT_TRUE(res.has_value());
        dir = dir.Norm();
        EXPECT_GE(res->Dot(normal), 0);
        EXPECT_LE(res->Dot(dir), 0);
    }
}



int main(int argc, char* argv[]){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

