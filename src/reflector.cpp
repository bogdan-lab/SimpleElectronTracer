#include "reflector.hpp"


std::optional<Vector> MirrorReflector::ReflectParticle(const Particle& pt,
                              const Vector& normal, std::mt19937& rnd_gen)const{
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    if(rnd(rnd_gen)>reflection_coefficient_){
        return std::nullopt; 	//particle died
    }
    double vel_proj = pt.GetDirection().Dot(normal);
    Vector new_vel = pt.GetDirection() - normal.Times(2.0*vel_proj);
    return new_vel.Norm();
}

std::optional<Vector> LambertianReflector::ReflectParticle(const Particle& pt,
                 const Vector& normal, std::mt19937& rnd_gen)const{
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    if(rnd(rnd_gen)>reflection_coefficient_){
        return std::nullopt; 	//particle died
    }
    return pt.GetRandomVel(normal, rnd_gen);
}
