
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


Particle::Particle(const Vector& given_p, const Vector& given_v){
    pos_ = given_p;
    V_ = given_v;
    vol_count_ = 0;
    surf_count_ = 0;
}

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
