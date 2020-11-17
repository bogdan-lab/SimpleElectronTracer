#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <string>
#include <memory>
#include <iostream>
#include <fstream>

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
    std::vector<Particle> stat_;
    std::vector<Vec3> contour_; 	//points which build the surface contour
    std::unique_ptr<Reflector> reflector_;
    SurfaceCoeficients coefs_;
    Vec3 normal_;
    FILE* output_file_;

    static SurfaceCoeficients CalcSurfaceCoefficients(
                                              const std::vector<Vec3> contour);
    void SaveSurfaceParticles() const;
public:

    Surface() = delete;
    Surface(std::vector<Vec3> g_contour,
            std::unique_ptr<Reflector> g_reflector, FILE* out_file,
            size_t dump_size);
    void WriteFileHeader() const;
    void SaveParticle(Particle&& pt);
    bool CheckIfPointOnSurface(const Vec3& point) const;
    std::vector<Vec3> TranslateContourIntoBasis(const ONBasis_3x3& basis) const;

    Vec3 GetPointOnSurface() const;
    const std::vector<Vec3>& GetContour() const ;
    const Vec3& GetNormal() const;
    bool IsSaveStat() const;
    const Reflector* GetReflector() const ;
    const SurfaceCoeficients& GetSurfaceCoefficients() const ;
    ~Surface();
};



#endif
