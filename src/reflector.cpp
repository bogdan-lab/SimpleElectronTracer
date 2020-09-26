#include "reflector.hpp"


std::pair<bool, Vector> MirrorReflector::ReflectParticle(const Particle& pt,
                              const Vector& normal, std::mt19937& rnd_gen)const{
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    if(rnd(rnd_gen)>reflection_coefficient_){
        return std::make_pair(false, Vector()); 	//particle died
    }
    double vel_proj = pt.GetDirection().Dot(normal);
    Vector new_vel = pt.GetDirection() - normal.Times(2.0*vel_proj);
    return std::make_pair(true, new_vel.Norm());
}

std::pair<bool, Vector> LambertianReflector::ReflectParticle(const Particle& pt,
                 const Vector& normal, std::mt19937& rnd_gen)const{
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    if(rnd(rnd_gen)>reflection_coefficient_){
        return std::make_pair(false, Vector()); 	//particle died
    }
    return std::make_pair(true, pt.GetRandomVel(normal, rnd_gen));
}
