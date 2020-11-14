#include <gtest/gtest.h>
#include "reflector.hpp"


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
        dir.Norm();
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
        dir.Norm();
        EXPECT_GE(res->Dot(normal), 0);
        EXPECT_LE(res->Dot(dir), 0);
    }
}
