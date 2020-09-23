
#include <random>
#include <utility>

#include "utils.hpp"
#include "surface.hpp"
#include "particle.hpp"




Particle::Particle(const Surface &s, default_random_engine &rnd_gen){
    pos_ = s.GetRandomPoint(rnd_gen);
    vol_count = 0;
    surf_count = 0;
    //rnd =  uniform_real_distribution<double>(0.0, 1.0);
    //V_ = this->GetRandVel(1, rnd_gen);     //directed along positive y
    V_ = {0.0, 1.0, 0.0}; //TODO decide about velocity generation!
}


Particle::Particle(const Point& given_p, const Velocity& given_v){
    pos_ = given_p;
    V_ = given_v;
    vol_count_ = 0;
    surf_count_ = 0;
}


std::pair<bool, Point> Particle::GetCrossPoint(const Surface& s) const {
    Point p_new(0.0, 0.0, 0.0);
    cross_flag = false;
    int coor_flag = s.GetCoorFlag();
    double coor_val = s.GetCoorVal();
    vector<double> x_bnd = s.GetXbnd();
    vector<double> y_bnd = s.GetYbnd();
    vector<double> z_bnd = s.GetZbnd();
    double t = 0.0;
    if (coor_flag==0){
        if((p.x<coor_val && V[0]>0.0) || (p.x>coor_val && V[0]<0.0)){
            //direction is correct
            p_new.x = coor_val;
            t = (coor_val - p.x)/V[0];
            p_new.y = p.y + V[1]*t;
            p_new.z = p.z + V[2]*t;
            if ((p_new.z>z_bnd[0] && p_new.z<z_bnd[1]) && (p_new.y>y_bnd[0] && p_new.y<y_bnd[1])){
                cross_flag = true;
            }
            else {
                cross_flag = false;
            }
        }
    }
    else if (coor_flag==1){
        if((p.y<coor_val && V[1]>0.0) || (p.y>coor_val && V[1]<0.0)){
            //direction is correct
            p_new.y = coor_val;
            t = (coor_val - p.y)/V[1];
            p_new.x = p.x + V[0]*t;
            p_new.z = p.z + V[2]*t;
            if ((p_new.z>z_bnd[0] && p_new.z<z_bnd[1]) && (p_new.x>x_bnd[0] && p_new.x<x_bnd[1])){
                cross_flag = true;
            }
            else {
                cross_flag = false;
            }
        }
    }
    else{
        if((p.z<coor_val && V[2]>0.0) || (p.z>coor_val && V[2]<0.0)){
            //direction is correct
            p_new.z = coor_val;
            t = (coor_val - p.z)/V[2];
            p_new.y = p.y + V[1]*t;
            p_new.x = p.x + V[0]*t;
            if ((p_new.x>x_bnd[0] && p_new.x<x_bnd[1]) && (p_new.y>y_bnd[0] && p_new.y<y_bnd[1])){
                cross_flag = true;
            }
            else {
                cross_flag = false;
            }
        }
    }
    return p_new;
}
