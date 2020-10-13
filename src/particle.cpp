#include <optional>
#include <random>
#include <utility>
#include <limits>

#include "particle.hpp"



Particle::Particle(const Vector& given_p, const Vector& given_v){
    pos_ = given_p;
    V_ = given_v.Norm();
    vol_count_ = 0;
    surf_count_ = 0;
}

Particle::Particle(const Vector &given_p, const Vector& direction,
                   std::mt19937& rnd_gen){
    pos_ = given_p;
    V_ = GetRandomVel(direction, rnd_gen);
    vol_count_ = 0;
    surf_count_ = 0;
}

const Vector& Particle::GetPosition() const{return pos_;}
const Vector& Particle::GetDirection() const {return V_;}
size_t Particle::GetVolCount() const {return vol_count_;}
size_t Particle::GetSurfCount() const {return surf_count_;}



std::optional<Vector> Particle::GetCrossPoint(const Surface& s) const {
    const Surface::SurfaceCoeficients& Sc = s.GetSurfaceCoefficients();
    if((Sc.A_*V_.GetX() + Sc.B_*V_.GetY() + Sc.C_*V_.GetZ()) == 0.0){
        //particle moves parallel to the surface
        return std::nullopt;
    }
    //Look at time needed to reach the surface
    double t = -1*(Sc.A_*pos_.GetX() + Sc.B_*pos_.GetY()
                                                 + Sc.C_*pos_.GetZ() + Sc.D_)
                        /(Sc.A_*V_.GetX() + Sc.B_*V_.GetY() + Sc.C_*V_.GetZ());
    if(t<=0){
        return std::nullopt;
    }
    //Here at least direction is correct --> check for boundaries
    Vector cross_point = {pos_.GetX() + V_.GetX()*t,
                          pos_.GetY() + V_.GetY()*t,
                          pos_.GetZ() + V_.GetZ()*t};
    cross_point = VerifyPointOnSurface(s, cross_point);
    //TODO think about how to overcome double precision problem adequately
    if (s.CheckIfPointOnSurface(cross_point)){
        return cross_point;
    }
    return std::nullopt;
}

double Particle::GetDistanceInGas(const Background& gas,
                                  std::mt19937& rnd_gen) const{
    if (gas.p_ == 0.0){
        return std::numeric_limits<double>::max();
    }
    double mfp = 1.38e-17*gas.T_/(gas.p_*gas.sigma_);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    return mfp*log(1.0/(1.0-rnd(rnd_gen)));
}


void Particle::MakeGasCollision(const double distance,
                                std::mt19937& rnd_gen){
    pos_ = pos_ + V_.Times(distance);
    vol_count_++;
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    double costheta = 2*rnd(rnd_gen)-1;
    double sintheta = sqrt(1-costheta*costheta);
    double phi = rnd(rnd_gen)*2*M_PI;
    V_ = Vector(sintheta*sin(phi),
                sintheta*cos(phi),
                costheta);
}

Vector Particle::GetRandomVel(const Vector& direction,
                                  std::mt19937& rnd_gen) const{
    std::uniform_real_distribution<double> rnd(0.0, 1.0);
    double cos_theta;
    double sin_theta;
    double phi;
    cos_theta = rnd(rnd_gen);
    sin_theta = sqrt(1-cos_theta*cos_theta);
    phi = rnd(rnd_gen)*2*M_PI;
    std::vector<Vector> coor_transition = GenerateONBasisByNewZ(direction);
    Vector res_vec = ApplyCoordinateTransition(coor_transition,
                     Vector(sin_theta*sin(phi), sin_theta*cos(phi), cos_theta));
    return res_vec;
}


bool Particle::MakeStep(std::vector<Surface>& walls,
                        const Background& gas, std::mt19937 &rnd_gen){
    double min_dist = GetDistanceInGas(gas, rnd_gen);
    size_t wall_id = 0;
    Vector point_on_surf;
    bool colide_in_gas_flag = true;
    for(size_t i=0; i<walls.size(); i++){
        auto cross_res = GetCrossPoint(walls[i]);
        if(cross_res && GetDistance(pos_, cross_res.value())<min_dist){
            min_dist = GetDistance(pos_, cross_res.value());
            wall_id = i;
            colide_in_gas_flag = false;
            point_on_surf = cross_res.value();
        }
    }
    if(colide_in_gas_flag){
        MakeGasCollision(min_dist, rnd_gen);
        return true;
    }
    //Here we collide with surface --> can die
    pos_ = point_on_surf;
    surf_count_++;
    auto surf_refl = walls[wall_id].GetReflector()->ReflectParticle(*this,
                                           walls[wall_id].GetNormal(), rnd_gen);
    if(surf_refl.first){
        V_ = surf_refl.second;
        return true;
    }
    //Here particle is dead --> save its position
    if (walls[wall_id].IsSaveStat()){
        walls[wall_id].SaveParticle(*this);
    }
    return false;
}





