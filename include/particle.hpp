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
    Vec3 pos_ = {};
    Vec3 V_ = {};
    size_t vol_count_ = {}; 	//number of volume collisions happened
    size_t surf_count_ = {};	//number of surface collisions happened
public:
    Particle() = default;
    Particle(const Vec3& given_p, const Vec3& given_v);
    Particle(const Vec3& given_p, const Vec3& direction,
                                                     std::mt19937& rnd_gen);
    std::optional<Vec3> GetCrossPoint(const Surface& s) const;
    double GetDistanceInGas(const Background& gas,
                            std::mt19937& rnd_gen) const;
    void MakeGasCollision(const double distance,
                          std::mt19937& rnd_gen);
    bool MakeStep(std::vector<Surface>& walls, const Background& gas,
                  std::mt19937& rnd_gen);
    Vec3 GetRandomVel(const Vec3& direction, std::mt19937& rnd_gen) const;

    const Vec3& GetPosition() const;
    const Vec3& GetDirection() const;
    size_t GetVolCount() const;
    size_t GetSurfCount() const;
};


#endif
