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

class Surface{
private:
    std::list<Particle> stat_;
    bool save_stat_;             //flag for saving statistics
    std::vector<Vector> contour_; 	//points which build the surface contour
    std::string surface_name_;
    std::unique_ptr<Reflector> reflector_;
    SurfaceCoeficients coefs_;
    Vector normal_;

    static SurfaceCoeficients CalcSurfaceCoefficients(const std::vector<Point> contour);
public:

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

    const std::vector<Vector>& GetContour() const ;
    bool IsSaveStat() const;
    void SaveParticle(const Particle& pt);
    //TODO some function about reflection???
};



#endif
