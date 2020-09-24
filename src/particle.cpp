
#include <random>
#include <utility>

#include "utils.hpp"
#include "surface.hpp"
#include "particle.hpp"




Particle::Particle(const Vector& given_p, const Vector& given_v){
    pos_ = given_p;
    V_ = given_v;
    vol_count_ = 0;
    surf_count_ = 0;
}

const Vector& Particle::GetPosition() const{return pos_;}
const Vector& Particle::GetDirection() const {return V_;}
size_t Particle::GetVolCount() const {return vol_count_;}
size_t Particle::GetSurfCount() const {return surf_count_;}



std::pair<bool, Vector> Particle::GetCrossPoint(const Surface& s) const {
    Surface::SurfaceCoeficients Sc = s.GetSurfaceCoefficients();
    //Look at time neede to reach the surface
    double t = -1*(Sc.A_*pos_.GetX() + Sc.B_*pos_.GetY()
                                                 + Sc.C_*pos_.GetZ() + Sc.D_)
                        /(Sc.A_*V_.GetX() + Sc.B_*V_.GetY() + Sc.C_*V_.GetZ());
    if(t<0){
        return std::make_pair(false, Vector(0.0, 0.0, 0.0));
    }
    //Here at least direction is correct --> check for boundaries
    Vector cross_point = {pos_.GetX() + V_.GetX()*t,
                          pos_.GetY() + V_.GetY()*t,
                          pos_.GetZ() + V_.GetZ()*t};
    //Building ray from cross point on the current plain
    const std::vector<Vector>& surface_contour =s.GetContour();
    Vector ray_direction = Vector(cross_point, surface_contour[0]) +
                            Vector(cross_point, surface_contour[1]);
    size_t contour_cross_count = 0;
    for(size_t i=0; i<surface_contour.size()-1; i++){
        //TODO FIND CROSSECTIONS WITH CONTOUR LINES
    }
    return std::make_pair(false, cross_point);
}

std::pair<bool, double> Particle::GetDistanceToSurface(const Surface &s) const {
    auto cross_result = GetCrossPoint(s);
    if (cross_result.first){
        return std::make_pair(true, Vector(pos_, cross_result.second).Length());
    }
    return std::make_pair(false, -1.0);
}


double Particle::GetDistanceInGas(const Background& gas,
                                  default_random_engine& rnd_gen) const{
    if (gas.p_ == 0.0){
        return -1.0;
    }
    double mfp = 1.38e-17*gas.T_/(gas.p_*gas.sigma_);
    uniform_real_distribution<double> rnd(0.0, 1.0);
    return mfp*log(1.0/(1.0-rnd(rnd_gen)));
}


void Particle::MakeGasCollision(const double distance,
                                default_random_engine& rnd_gen){
    pos_ += distance*V_;
    vol_count_++;
    uniform_real_distribution<double> rnd(0.0, 1.0);
    double costheta = 2*rnd(rnd_gen)-1;
    double sintheta = sqrt(1-costheta*costheta);
    double phi = rnd(rdn_gen)*2*M_PI;
    V_ = Vector(sintheta*sin(phi),
                sintheta*cos(phi),
                costheta);
}

Vector Particle::GetRandomVel(const Vector &direction,
                              default_random_engine &rnd_gen){
    //TODO generate random velocity in hemisphere according to the given direction
    return false;
}


bool Particle::MakeStep(const std::vector<Surface>& walls,
                        const Background& gas, default_random_engine &rnd_gen){
    double min_dist = GetDistanceInGas(gas, rnd_gen);
    size_t wall_id = 0;
    Vector point_on_surf(0.0, 0.0, 0.0);
    bool colide_in_gas_flag = true;
    for(size_t i=0; i<walls.size(); i++){
        auto cross_res = GetCrossPoint(walls[i]);
        if(cross_res.first && GetDistance(pos_, cross_res.second)<min_dist){
            min_dist = GetDistance(pos_, cross_res.second);
            wall_id = i;
            colide_in_gas_flag = false;
            point_on_surf = cross_res.second;
        }
    }
    if(!colide_in_gas_flag){
        MakeGasCollision(min_dist, rnd_gen);
        return true;
    }
    //Here we collide with surface --> can die
    pos_ = point_on_surf;
    surf_count_++;
    //TODO change velocity if i live or die
}





