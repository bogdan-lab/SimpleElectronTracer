#ifndef REFLECTOR_HEADER
#define REFLECTOR_HEADE

#include <random>

#include "particle.hpp"
#include "surface.hpp"

class Reflector{
public:
    virtual bool ReflectParticle(Particle& pt, Surface& s, default_random_engine& rnd_gen) const = 0;
    virtual ~Reflector() = default;
};

class UniformReflector : public Reflector {
private:
    double reflection_coefficient_;
public:
    UniformReflector(const double val): reflection_coefficient_(val) {}
    bool ReflectParticle(Particle &pt, Surface &s, default_random_engine &rnd_gen) const override;
};

class LambertianReflector : public Reflector {
private:
    double reflection_coefficient_;
public:
    LambertianReflector(const double val): reflection_coefficient_(val) {}
    bool ReflectParticle(Particle &pt, Surface &s, default_random_engine &rnd_gen) const override;
}


#endif //REFLECTOR_HEADER
