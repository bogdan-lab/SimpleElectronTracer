#include <gtest/gtest.h>
#include "surface.hpp"

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
