#ifndef SURFACE_HEADER
#define SURFACE_HEADER

#include <string>
#include <memory>

#include "particle.hpp"
#include "reflector.hpp"

class Surface{
private:
    std::vector<Particle> stat_;
    bool save_stat_;             //flag for saving statistics
    std::vector<Vector> contour_; 	//points which build the surface contour
    std::string surface_name_;
    std::unique_ptr<Reflector> reflector_;
    SurfaceCoeficients coefs_;

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
    const SurfaceCoeficients& GetSurfaceCoefficients() const ;
    const std::vector<Vector>& GetContour() const ;


    Surface(std::string& line);
    void SaveInfo(const Particle& pt, const vector<Point>& pt_tragectory);
    void WriteInfo(const size_t number_of_simulated_particles, const double pressure);
    Point GetRandomPoint(default_random_engine& rnd_gen) const;
    void PrintSurface();

    int GetCoorFlag() const;
    double GetCoorVal() const;
    vector<double> GetXbnd() const ;
    vector<double> GetYbnd() const ;
    vector<double> GetZbnd() const ;
    double GetRefl() const ;
    bool GetSaveStatFlag();
};



#endif
