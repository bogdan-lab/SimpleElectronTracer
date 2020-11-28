#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <string>
#include <memory>
#include <iostream>
#include <fstream>

#include "particle.hpp"
#include "reflector.hpp"
#include "math.hpp"

class Reflector;
class Particle;

class Surface{
public:
    struct SurfaceCoeficients{
        //Ax + By + Cz + D = 0
        double A_;
        double B_;
        double C_;
        double D_;
    };

private:
    std::vector<Vec3> contour_; 	//points which build the surface contour
    std::unique_ptr<Reflector> reflector_;
    SurfaceCoeficients coefs_;
    std::vector<double> tri_areas_;
    Vec3 mass_center_;
    ONBasis_3x3 surf_basis_;
    std::vector<Vec3> basis_contour_;
    bool save_stat_flag_;

    int GetQuarter(const Vec3& point, const Vec3& node) const;
    double GetOrientationWinding(const Vec3& point, const Vec3& prev_node,
                                 const Vec3& next_node) const;
    int CalcWindChange(const Vec3& prev_node, const Vec3& next_node,
                            const Vec3& point) const ;

public:

    Surface(std::vector<Vec3>&& g_contour, std::unique_ptr<Reflector>&& g_reflector,
            bool g_save_stat_flag);
    bool CheckIfPointOnSurface(const Vec3& point) const;
    std::optional<Vec3> GetCrossPoint(const Vec3& position,
                                      const Vec3& direction) const;
    void VerifyPointInVolume(const Vec3& start, Vec3 &end) const;

    Vec3 GetRandomPointInContour(std::mt19937& rng) const;
    const Vec3& GetMassCenter() const;
    const std::vector<Vec3>& GetContour() const ;
    const Vec3& GetNormal() const;
    bool IsSaveStat() const;
    const Reflector* GetReflector() const ;
    const SurfaceCoeficients& GetSurfaceCoefficients() const ;

    static std::vector<double> CalcTriangleAreas(const std::vector<Vec3>& contour);
    static Vec3 CalcCenterOfMass(const std::vector<Vec3>& contour);
    static std::vector<Vec3> TranslateContourIntoBasis(const ONBasis_3x3& basis,
                                              const std::vector<Vec3>& contour);
    static SurfaceCoeficients CalcSurfaceCoefficients(
                                                const std::vector<Vec3> contour);
};



#endif
