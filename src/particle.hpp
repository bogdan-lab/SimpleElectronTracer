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
    Particle(const Vector& given_p, const Vector& given_v);
    std::pair<bool, Vector> GetCrossPoint(const Surface& s) const;
    std::pair<bool, double> GetDistanceToSurface(const Surface& s) const;
    double GetDistanceInGas(const Background& gas,
                            default_random_engine& rnd_gen) const;
    void MakeGasCollision(const double distance,
                          default_random_engine& rnd_gen);
    bool MakeStep(const std::vector<Surface>& walls, const Background& gas,
                  default_random_engine& rnd_gen);
    Vector GetRandomVel(const Vector& direction,
                        default_random_engine& rnd_gen);

    const Vector& GetPosition() const;
    const Vector& GetDirection() const;
    size_t GetVolCount() const;
    size_t GetSurfCount() const;
};


#endif
