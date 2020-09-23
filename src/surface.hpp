#ifndef SURFACE_HEADER
#define SURFACE_HEADER

#include <string>
#include <memory>

#include "particle.hpp"
#include "reflector.hpp"

class Surface{
public:
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
private:
    std::vector<Particle> stat_;
    bool save_stat_;             //flag for saving statistics
    std::vector<Point> contour_; 	//points which build the surface contour
    std::string surface_name_;
    std::unique_ptr<Reflector> reflector_;
};



#endif
