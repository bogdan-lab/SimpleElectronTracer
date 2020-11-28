#include <gtest/gtest.h>
#include "surface.hpp"

TEST(SurfaceTests, GenerationTest){
    std::ofstream file;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(std::move(contour), std::make_unique<MirrorReflector>(0.0),
              std::move(file));
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


TEST(SurfaceTests, CheckIfPointOnSurface){
    std::ofstream file;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    Surface s(std::move(contour), std::make_unique<MirrorReflector>(0.0),
              std::move(file));
    Vec3 point(1.0, 0.5, 0.7);
    EXPECT_TRUE(s.CheckIfPointOnSurface(point));
    Vec3 point1(1.0+0.1, 0.5, 0.7);
    EXPECT_TRUE(s.CheckIfPointOnSurface(point1));
    Vec3 point2(1.0, 12.0, 0.5);
    EXPECT_FALSE(s.CheckIfPointOnSurface(point2));
}


TEST(SurfaceTests, VerifyPointOnSurfaceTest){
    std::ofstream file;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    std::unique_ptr<Surface> s = std::make_unique<Surface>(std::move(contour),
           std::make_unique<MirrorReflector>(0.0), std::move(file));
    Vec3 end(1+2e-6, 0.5, 0.4);
    Vec3 start(0.5, 0.5, 0.4);
    s->VerifyPointInVolume(start, end);
    EXPECT_NEAR(end.GetX(), 1.0, 1e-15);
    EXPECT_LE(end.GetX(), 1.0);
    EXPECT_EQ(end.GetY(), 0.5);
    EXPECT_EQ(end.GetZ(), 0.4);

    Vec3 start2(0.5, 0.3, 0.0);
    Vec3 end2(1+1, 0.3, 1.5);
    s->VerifyPointInVolume(start2, end2);
    EXPECT_NEAR(end2.GetX(), 1.0, 1e-15);
    EXPECT_LE(end2.GetX(), 1.0);
    EXPECT_EQ(end2.GetY(), 0.3);
    EXPECT_NEAR(end2.GetZ(), 0.5, 1e-15);
}

TEST(SurfaceTests, TriangleAreas){
    std::vector<Vec3> contour = {Vec3(1.0, 0.0, 0.0),
                                 Vec3(1.0, 0.0, 1.0),
                                 Vec3(1.0, 1.0, 1.0),
                                 Vec3(1.0, 2.0, 0.0)};
    std::vector<double> areas = Surface::CalcTriangleAreas(contour);
    EXPECT_EQ(areas.size(), 2);
    EXPECT_EQ(areas[0], 0.5);
    EXPECT_EQ(areas[1], 1.0);
}


TEST(SurfaceTests, RandomPointGeneration){
    std::mt19937 rng(42u);
    std::ofstream file;
    std::unique_ptr<char[]> buff;
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 2.0, 0.0)};
    std::unique_ptr<Surface> s = std::make_unique<Surface>(std::move(contour),
           std::make_unique<MirrorReflector>(0.0), std::move(file));
    for(size_t i=0; i<100; i++){
        Vec3 point = s->GetRandomPointInContour(rng);
        EXPECT_TRUE(s->CheckIfPointOnSurface(point));
    }
}
