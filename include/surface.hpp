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
    struct Boundary{
        double min_;
        double max_;
    };

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
    //TODO Redo into fully plygonal representaiotn
    Boundary x_bnd_;
    Boundary y_bnd_;
    Boundary z_bnd_;

    static Boundary GetBoundary(const std::vector<Vec3>& ctr, const char axis);
    static SurfaceCoeficients CalcSurfaceCoefficients(
                                              const std::vector<Vec3> contour);
public:


    Surface(std::string g_name, std::vector<Vec3> g_contour,
            std::unique_ptr<Reflector> g_reflector, bool save_flag);
    void SaveSurfaceParticles() const;
    void SaveParticle(const Particle& pt);
    bool CheckIfPointOnSurface(const Vec3& point) const;

    const std::vector<Vec3>& GetContour() const ;
    const Vec3& GetNormal() const;
    bool IsSaveStat() const;
    const Reflector* GetReflector() const ;
    const std::string& GetName() const ;
    const SurfaceCoeficients& GetSurfaceCoefficients() const ;
    const Boundary& GetXBnd() const ;
    const Boundary& GetYBnd() const ;
    const Boundary& GetZBnd() const ;
};



#endif
