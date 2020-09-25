#ifndef REFLECTOR_HEADER
#define REFLECTOR_HEADE

#include <random>

#include "particle.hpp"
#include "surface.hpp"

class Surface;

class Reflector{
public:
    //TODO Make it return second optional for the case when particle dies!
    virtual std::pair<bool, Vector> ReflectParticle(const Particle& pt,
                           const Vector& normal, std::mt19937& rnd_gen) const = 0;
    virtual ~Reflector() = default;
};

class MirrorReflector : public Reflector {
private:
    double reflection_coefficient_;
public:
    MirrorReflector(const double val): reflection_coefficient_(val) {}
    std::pair<bool, Vector> ReflectParticle(const Particle &pt,
                  const Vector& normal, std::mt19937& rnd_gen) const override;
};

class LambertianReflector : public Reflector {
private:
    double reflection_coefficient_;
public:
    LambertianReflector(const double val): reflection_coefficient_(val) {}
    std::pair<bool, Vector> ReflectParticle(const Particle &pt,
                 const Vector& normal, std::mt19937& rnd_gen) const override;
}


#endif //REFLECTOR_HEADER
