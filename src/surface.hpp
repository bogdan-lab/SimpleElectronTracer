#ifndef SURFACE_HEADER
#define SURFACE_HEADER

#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <list>

#include "particle.hpp"
#include "reflector.hpp"
#include "utils.hpp"

class Reflector;

class Surface{
private:
    std::list<Particle> stat_;
    bool save_stat_;             //flag for saving statistics
    std::vector<Vector> contour_; 	//points which build the surface contour
    std::string surface_name_;
    std::unique_ptr<Reflector> reflector_;
    SurfaceCoeficients coefs_;
    Vector normal_;
    //TODO Redo into fully plygonal representaiotn
    Boundary x_bnd_;
    Boundary y_bnd_;
    Boundary z_bnd_;

    static Boundary GetBoundary(const std::vector<Vector>& ctr, const char axis);
    static SurfaceCoeficients CalcSurfaceCoefficients(
                                              const std::vector<Point> contour);
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

    Surface(std::string g_name, std::vector<Vector> g_contour,
            std::unique_ptr<Reflector> g_reflector, bool save_flag);
    void SaveSurfaceParticles(std::ofstream& out) const;
    void SaveParticle(const Particle& pt);
    bool CheckIfPointOnSurface(const Vector& point) const;

    const std::vector<Vector>& GetContour() const ;
    const Vector& GetNormal() const;
    bool IsSaveStat() const;
    const Reflector* GetReflector() const ;
    const std::string& GetName() const ;
};



#endif
