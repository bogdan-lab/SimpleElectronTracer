﻿#include <gtest/gtest.h>
#include "surface.hpp"

TEST(SurfaceTests, GenerationTest){
    std::ofstream file;
    std::unique_ptr<char[]> buff;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(contour, std::make_unique<MirrorReflector>(0.0), std::move(file),
               std::move(buff));
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
    std::ofstream file;
    std::unique_ptr<char[]> buff;
    std::vector<Vec3> contour {Vec3(0.0, 0.0, 0.0),
                               Vec3(0.5, 0.5, 0.0),
                               Vec3(0.5, 0.5, 1.0),
                               Vec3(0.0, 0.0, 1.0)};
    Surface s(contour, std::make_unique<MirrorReflector>(0.0), std::move(file),
               std::move(buff));
    Vec3 p = s.GetPointOnSurface();
    Surface::SurfaceCoeficients Sc=s.GetSurfaceCoefficients();
    EXPECT_EQ(p.GetX()*Sc.A_ + p.GetY()*Sc.B_ + p.GetZ()*Sc.C_ + Sc.D_, 0.0);
}

TEST(SurfaceTests, CheckIfPointOnSurface){
    std::ofstream file;
    std::unique_ptr<char[]> buff;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(contour, std::make_unique<MirrorReflector>(0.0), std::move(file),
               std::move(buff));
    Vec3 point(1.0, 0.5, 0.7);
    EXPECT_TRUE(s.CheckIfPointOnSurface(point));
    Vec3 point1(1.0+0.1, 0.5, 0.7);
    EXPECT_TRUE(s.CheckIfPointOnSurface(point1));
    Vec3 point2(1.0, 12.0, 0.5);
    EXPECT_FALSE(s.CheckIfPointOnSurface(point2));
}


TEST(SurfaceTests, VerifyPointOnSurfaceTest){
    std::ofstream file;
    std::unique_ptr<char[]> buff;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    std::unique_ptr<Surface> s = std::make_unique<Surface>(contour,
           std::make_unique<MirrorReflector>(0.0), std::move(file),
            std::move(buff));
    Vec3 end(1+2e-6, 0.5, 0.4);
    Vec3 start(0.5, 0.5, 0.4);
    //s->VerifyPointInVolume(start, end);
    //EXPECT_NEAR(end.GetX(), 1.0, 1e-15);
    //EXPECT_LE(end.GetX(), 1.0);
    //EXPECT_EQ(end.GetY(), 0.5);
    //EXPECT_EQ(end.GetZ(), 0.4);

    Vec3 start2(0.5, 0.3, 0.0);
    Vec3 end2(1+1, 0.3, 1.5);
    s->VerifyPointInVolume(start2, end2);
    EXPECT_NEAR(end2.GetX(), 1.0, 1e-15);
    EXPECT_LE(end2.GetX(), 1.0);
    EXPECT_EQ(end2.GetY(), 0.3);
    EXPECT_NEAR(end2.GetZ(), 0.5, 1e-15);
}

