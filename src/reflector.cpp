#include "reflector.hpp"


bool UniformReflector::ReflectParticle(Particle &pt, Surface &s, default_random_engine &rnd_gen) const{
    return false;
}

bool LambertianReflector::ReflectParticle(Particle &pt, Surface &s, default_random_engine &rnd_gen) const {
    return false;
}
