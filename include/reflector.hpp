#ifndef REFLECTOR_HPP
#define REFLECTOR_HPP

#include <random>

#include "particle.hpp"
#include "surface.hpp"
#include "utils.hpp"

class Surface;
class Particle;

class Reflector{
public:
    virtual std::optional<Vec3> ReflectParticle(const Particle& pt,
                           const Vec3& normal, std::mt19937& rnd_gen) const = 0;
    virtual ~Reflector() = default;
};

class MirrorReflector : public Reflector {
private:
    double reflection_coefficient_;
public:
    MirrorReflector(const double val): reflection_coefficient_(val) {}
    std::optional<Vec3> ReflectParticle(const Particle &pt,
                  const Vec3& normal, std::mt19937& rnd_gen) const override;
};

class LambertianReflector : public Reflector {
private:
    double reflection_coefficient_;
public:
    LambertianReflector(const double val): reflection_coefficient_(val) {}
    std::optional<Vec3> ReflectParticle(const Particle &pt,
                 const Vec3& normal, std::mt19937& rnd_gen) const override;
};


#endif //REFLECTOR_HEADER
