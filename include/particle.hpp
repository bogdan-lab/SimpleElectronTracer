#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>
#include <random>
#include <utility>
#include <optional>

#include "surface.hpp"
#include "utils.hpp"

class Surface;

class Particle{
private:
    Vector pos_;
    Vector V_;
    size_t vol_count_; 	//number of volume collisions happened
    size_t surf_count_;	//number of surface collisions happened
public:
    Particle(const Vector& given_p, const Vector& given_v);
    Particle(const Vector& given_p, const Vector& direction,
                                                     std::mt19937& rnd_gen);
    std::optional<Vector> GetCrossPoint(const Surface& s) const;
    std::pair<bool, double> GetDistanceToSurface(const Surface& s) const;
    double GetDistanceInGas(const Background& gas,
                            std::mt19937& rnd_gen) const;
    void MakeGasCollision(const double distance,
                          std::mt19937& rnd_gen);
    bool MakeStep(std::vector<Surface>& walls, const Background& gas,
                  std::mt19937& rnd_gen);
    Vector GetRandomVel(const Vector& direction, std::mt19937& rnd_gen) const;

    const Vector& GetPosition() const;
    const Vector& GetDirection() const;
    size_t GetVolCount() const;
    size_t GetSurfCount() const;
};


#endif
