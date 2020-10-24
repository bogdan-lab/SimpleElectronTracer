#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include "utils.hpp"
#include "particle.hpp"
#include "reflector.hpp"
#include "surface.hpp"

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


TEST(UtilsTests, VerifyPointOnSurfaceTest){
    std::string name = "test_surface";
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(name, contour, std::make_unique<MirrorReflector>(0.0), false);
    Vec3 point(1+2e-6, 0.5, 0.4);
    VerifyPointInVolume(s, point);
    EXPECT_NEAR(point.GetX(), 1.0, 1e-15);
    EXPECT_EQ(point.GetY(), 0.5);
    EXPECT_EQ(point.GetZ(), 0.4);
    Vec3 point2(0.9, 0.3, 0.4);
    EXPECT_EQ(point2.GetX(), 0.9);
    EXPECT_EQ(point2.GetY(), 0.3);
    EXPECT_EQ(point2.GetZ(), 0.4);
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

TEST(SurfaceTests, GenerationTest){
    std::string name = "test_surface";
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(name, contour, std::make_unique<MirrorReflector>(0.0), false);
    EXPECT_EQ(s.GetName(), name);
    EXPECT_NEAR(s.GetNormal().Length(), 1.0, 1e-15);
    EXPECT_EQ(s.GetNormal().GetX(), -1.0);
    EXPECT_EQ(s.GetNormal().GetY(), 0.0);
    EXPECT_EQ(s.GetNormal().GetZ(), 0.0);
    Surface::SurfaceCoeficients Sc = s.GetSurfaceCoefficients();
    EXPECT_EQ(Sc.A_, -1.0);
    EXPECT_EQ(Sc.B_, 0.0);
    EXPECT_EQ(Sc.C_, 0.0);
    EXPECT_EQ(Sc.D_, 1.0);
    EXPECT_FALSE(s.IsSaveStat());
}

TEST(SurfaceTests, GetPointOnSurfaceTest){
    std::string name = "test_surface";
    std::vector<Vec3> contour {Vec3(0.0, 0.0, 0.0),
                               Vec3(0.5, 0.5, 0.0),
                               Vec3(0.5, 0.5, 1.0),
                               Vec3(0.0, 0.0, 1.0)};
    Surface s(name, contour, std::make_unique<MirrorReflector>(0.0), false);
    Vec3 p = s.GetPointOnSurface();
    Surface::SurfaceCoeficients Sc=s.GetSurfaceCoefficients();
    EXPECT_EQ(p.GetX()*Sc.A_ + p.GetY()*Sc.B_ + p.GetZ()*Sc.C_ + Sc.D_, 0.0);
}

TEST(SurfaceTests, CheckIfPointOnSurface){
    std::string name = "test_surface";
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(name, contour, std::make_unique<MirrorReflector>(0.0), false);
    Vec3 point(1.0, 0.5, 0.7);
    EXPECT_TRUE(s.CheckIfPointOnSurface(point));
    Vec3 point1(1.0+0.1, 0.5, 0.7);
    EXPECT_TRUE(s.CheckIfPointOnSurface(point1));
    Vec3 point2(1.0, 12.0, 0.5);
    EXPECT_FALSE(s.CheckIfPointOnSurface(point2));
}


int main(int argc, char* argv[]){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

