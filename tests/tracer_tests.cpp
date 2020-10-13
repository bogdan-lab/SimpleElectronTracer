#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include "utils.hpp"
#include "particle.hpp"

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


TEST(UtilsTests, CoordinateTransitionTest){
    std::vector<Vector> basis(3);
    basis[0] = Vector(1.0, 0.0, 0.0);
    basis[1] = Vector(0.0, 1.0, 0.0);
    basis[2] = Vector(0.0, 0.0, -1.0);
    Vector tst(0.1, 0.2, 0.3);
    Vector res = ApplyCoordinateTransition(basis, tst);
    EXPECT_EQ(res.GetX(), tst.GetX());
    EXPECT_EQ(res.GetY(), tst.GetY());
    EXPECT_EQ(res.GetZ(), -tst.GetZ());
}


TEST(UtilsTests, ONBGenerationTest){
    Vector new_z(0.1, 0.2, 0.3);
    std::vector<Vector> basis = GenerateONBasisByNewZ(new_z);
    EXPECT_NEAR(1.0, basis[0].Length(), 1e-15);
    EXPECT_NEAR(1.0, basis[1].Length(), 1e-15);
    EXPECT_NEAR(1.0, basis[2].Length(), 1e-15);

    EXPECT_NEAR(0.0, basis[0].Dot(basis[1]), 1e-15);
    EXPECT_NEAR(0.0, basis[1].Dot(basis[2]), 1e-15);
    EXPECT_NEAR(0.0, basis[0].Dot(basis[2]), 1e-15);

    Vector z_orig(0.0, 0.0, 1.0);
    Vector z_trans = ApplyCoordinateTransition(basis, z_orig);
    new_z = new_z.Norm();
    EXPECT_NEAR(new_z.GetX(), z_trans.GetX(), 1e-15);
    EXPECT_NEAR(new_z.GetY(), z_trans.GetY(), 1e-15);
    EXPECT_NEAR(new_z.GetZ(), z_trans.GetZ(), 1e-15);
}


TEST(ParticleTests, ParticleGenerationtest){
    {
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
    {
        Particle pt(Vector(1.0, 2.0, 3.0), Vector(4.0, 5.0, 6.0));
        EXPECT_EQ(pt.GetPosition().GetX(), 1.0);
        EXPECT_EQ(pt.GetPosition().GetY(), 2.0);
        EXPECT_EQ(pt.GetPosition().GetZ(), 3.0);
        EXPECT_NEAR(pt.GetDirection().GetX(), 4.0/sqrt(77.0), 1e-15);
        EXPECT_NEAR(pt.GetDirection().GetY(), 5.0/sqrt(77.0), 1e-15);
        EXPECT_NEAR(pt.GetDirection().GetZ(), 6.0/sqrt(77.0), 1e-15);
        EXPECT_EQ(pt.GetSurfCount(), 0);
        EXPECT_EQ(pt.GetVolCount(), 0);
    }
    {
        std::mt19937 rnd_gen(42);
        Vector pos(-1.0, 15.0, 48.0);
        Vector dir(-2.0, 5.0, 10.0);
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
}

TEST(ParticleTests, GetRandomVelTest){
    std::mt19937 rnd_gen(42);
    Vector dir(1.0, -5.0, 8.0);
    Particle pt;
    for(size_t i=0; i<100; i++){
        Vector res = pt.GetRandomVel(dir, rnd_gen);
        EXPECT_GE(dir.Dot(res), 0);
        EXPECT_NEAR(res.Length(), 1.0, 1e-15);
    }
}

TEST(ParticleTests, MakeGasCollisionTest){
    double distance = 15;
    Vector dir(4.0, 5.0, 6.0);
    Vector pos(1.0, 2.0, 3.0);
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

int main(int argc, char* argv[]){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

