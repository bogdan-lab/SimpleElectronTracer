#include <gtest/gtest.h>

#include "utils.hpp"
#include "surface.hpp"

TEST(UtilsTests, VerifyPointOnSurfaceTest){
    std::vector<Vec3> contour {Vec3(1.0, 0.0, 0.0),
                               Vec3(1.0, 0.0, 1.0),
                               Vec3(1.0, 1.0, 1.0),
                               Vec3(1.0, 1.0, 0.0)};
    std::unique_ptr<Surface> s = std::make_unique<Surface>(contour, std::make_unique<MirrorReflector>(0.0), nullptr,5);
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
