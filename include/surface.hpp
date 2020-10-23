#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <list>

#include "particle.hpp"
#include "reflector.hpp"
#include "utils.hpp"

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
    std::list<Particle> stat_;
    bool save_stat_;             //flag for saving statistics
    std::vector<Vec3> contour_; 	//points which build the surface contour
    std::string surface_name_;
    std::unique_ptr<Reflector> reflector_;
    SurfaceCoeficients coefs_;
    Vec3 normal_;

    static SurfaceCoeficients CalcSurfaceCoefficients(
                                              const std::vector<Vec3> contour);
public:

    Surface() = delete;
    Surface(std::string g_name, std::vector<Vec3> g_contour,
            std::unique_ptr<Reflector> g_reflector, bool save_flag);
    void SaveSurfaceParticles() const;
    void SaveParticle(const Particle& pt);
    std::optional<Vec3> CheckIfPointOnSurface(const Vec3& point) const;
    std::vector<Vec3> TranslateContourIntoBasis(const ONBasis_3x3& basis) const;

    Vec3 GetPointOnSurface() const;
    const std::vector<Vec3>& GetContour() const ;
    const Vec3& GetNormal() const;
    bool IsSaveStat() const;
    const Reflector* GetReflector() const ;
    const std::string& GetName() const ;
    const SurfaceCoeficients& GetSurfaceCoefficients() const ;
};



#endif
