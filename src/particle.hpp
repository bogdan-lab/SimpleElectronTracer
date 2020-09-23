#ifndef PARTICLE_HEADER
#define PARTICLE_HAEDER

#include <vector>
#include <random>
#include <utility>

#include "utils.hpp"
#include "surface.hpp"

class Particle{
private:
    Vector pos_;
    Vector V_;
    size_t vol_count_; 	//number of volume collisions happened
    size_t surf_count_;	//number of surface collisions happened
public:
    Particle(const Surface& s, default_random_engine& rnd_gen);
    Particle(const Vector& given_p, const Vector& given_v);
    std::pair<bool, Vector> GetCrossPoint(const Surface& s) const;



    double GetDistanceToSurface(const Surface& s);
    int GetReflectionSurfaceID(const std::vector<Surface>& walls);
    bool ReflectSurface(Surface& s, default_random_engine& rnd_gen);
    void MakeGasCollision(double pt_dist, default_random_engine& rnd_gen);
    Point GetPosition() const;
    double GetDistanceInGas(const double pressure, default_random_engine& rnd_gen) const;
    vector<double>GetRandVel(int direction, default_random_engine& rnd_gen) const;

};


#endif
