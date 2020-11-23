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
    std::vector<Vec3> contour_; 	//points which build the surface contour
    std::unique_ptr<Reflector> reflector_;
    std::unique_ptr<char[]> io_buffer_;
    std::ofstream output_file_;
    SurfaceCoeficients coefs_;
    Vec3 normal_;

    static SurfaceCoeficients CalcSurfaceCoefficients(
                                              const std::vector<Vec3> contour);
public:

    Surface(std::vector<Vec3> g_contour,
            std::unique_ptr<Reflector> g_reflector, std::ofstream&& out_file,
             std::unique_ptr<char[]>&& buff);
    void WriteFileHeader();
    void SaveParticle(Particle&& pt);
    bool CheckIfPointOnSurface(const Vec3& point) const;
    std::vector<Vec3> TranslateContourIntoBasis(const ONBasis_3x3& basis) const;
    std::optional<Vec3> GetCrossPoint(const Vec3& position,
                                      const Vec3& direction) const;
    void VerifyPointInVolume(const Vec3& start, Vec3 &end,
                             const double epsilon=1e-15) const;

    Vec3 GetPointOnSurface() const;
    const std::vector<Vec3>& GetContour() const ;
    const Vec3& GetNormal() const;
    bool IsSaveStat() const;
    const Reflector* GetReflector() const ;
    const SurfaceCoeficients& GetSurfaceCoefficients() const ;
};



#endif
