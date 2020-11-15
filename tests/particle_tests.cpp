#include <gtest/gtest.h>
#include "particle.hpp"


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
    dir.Norm();
    EXPECT_NEAR(pt.GetPosition().GetX(), pos.GetX() + distance*dir.GetX(), 1e-15);
    EXPECT_NEAR(pt.GetPosition().GetY(), pos.GetY() + distance*dir.GetY(), 1e-15);
    EXPECT_NEAR(pt.GetPosition().GetZ(), pos.GetZ() + distance*dir.GetZ(), 1e-15);
    EXPECT_NE(dir.GetX(), pt.GetDirection().GetX());
    EXPECT_NE(dir.GetY(), pt.GetDirection().GetY());
    EXPECT_NE(dir.GetZ(), pt.GetDirection().GetZ());
}


TEST(ParticleTests, GeneratorTest){
    auto rand_pt_gen = Particle::GetGenerator(true);
    auto stat_pt_gen = Particle::GetGenerator(false);
    std::mt19937 rnd_gen;
    rnd_gen.seed(42);
    Vec3 start_point(0, 1, 2);
    Vec3 direction(5,6,7);
    Particle stat_pt = stat_pt_gen(start_point, direction, rnd_gen);
    EXPECT_TRUE(stat_pt.GetDirection()== direction.Norm());
    Vec3 direction2(5,6,7);
    Particle rnd_pt = rand_pt_gen(start_point, direction2, rnd_gen);
    EXPECT_NE(rnd_pt.GetDirection(), direction2.Norm());
}
