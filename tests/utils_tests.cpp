#include <gtest/gtest.h>

#include "utils.hpp"
#include "surface.hpp"


TEST(UtilsTests, ChekCorrectGeometryTest){
    std::vector<std::unique_ptr<Surface>> walls;
    std::ofstream out_f_1;
    std::unique_ptr<char[]> no_ptr_1;
    std::unique_ptr<Reflector> ref_1 = std::make_unique<LambertianReflector>(0.1);
    std::vector<Vec3> contour_1 ={Vec3{1.0, 0.0, 0.0},
                                Vec3{1.0, 0.0, 1.0},
                                Vec3{1.0, 1.0, 1.0},
                                Vec3{1.0, 1.0, 0.0}};
    walls.push_back(std::make_unique<Surface>(std::move(contour_1), std::move(ref_1),
                                     std::move(out_f_1), std::move(no_ptr_1)));

    std::ofstream out_f_2;
    std::unique_ptr<char[]> no_ptr_2;
    std::unique_ptr<Reflector> ref_2 = std::make_unique<LambertianReflector>(0.1);
    std::vector<Vec3> contour_2 = {Vec3{0.0, 0.0, 0.0},
                                Vec3{0.0, 1.0, 0.0},
                                Vec3{0.0, 1.0, 1.0},
                                Vec3{0.0, 0.0, 1.0}};
    walls.push_back(std::make_unique<Surface>(std::move(contour_2), std::move(ref_2),
                                     std::move(out_f_2), std::move(no_ptr_2)));
    EXPECT_TRUE(check_surface_orientations(walls));

}

TEST(UtilsTests, ChekWrongGeometryTest){
    std::vector<std::unique_ptr<Surface>> walls;
    std::ofstream out_f_1;
    std::unique_ptr<char[]> no_ptr_1;
    std::unique_ptr<Reflector> ref_1 = std::make_unique<LambertianReflector>(0.1);
    std::vector<Vec3> contour_1 = {Vec3{1.0, 0.0, 0.0},
                                Vec3{1.0, 0.0, 1.0},
                                Vec3{1.0, 1.0, 1.0},
                                Vec3{1.0, 1.0, 0.0}};
    walls.push_back(std::make_unique<Surface>(std::move(contour_1), std::move(ref_1),
                                     std::move(out_f_1), std::move(no_ptr_1)));

    std::ofstream out_f_2;
    std::unique_ptr<char[]> no_ptr_2;
    std::unique_ptr<Reflector> ref_2 = std::make_unique<LambertianReflector>(0.1);
    std::vector<Vec3> contour_2 = {Vec3{0.0, 0.0, 0.0},
                                Vec3{0.0, 0.0, 1.0},
                                Vec3{0.0, 1.0, 1.0},
                                Vec3{0.0, 1.0, 0.0}};
    walls.push_back(std::make_unique<Surface>(std::move(contour_2), std::move(ref_2),
                                     std::move(out_f_2), std::move(no_ptr_2)));
    EXPECT_FALSE(check_surface_orientations(walls));

}
