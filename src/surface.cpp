
#include <utility>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>

#include "surface.hpp"

Surface::Surface(std::string g_name, std::vector<Vector> g_contour,
        std::unique_ptr<Reflector> g_reflector, bool save_flag):
    stat_(), save_stat_(save_flag), contour_(std::move(g_contour)),
    surface_name_(std::move(g_name)), reflector_(std::move(g_reflector))
{
    if(contour_.size()!=4){
        std::cerr << "CURRENT VERSION WORKS ONLY WITH RECTANGLE POLYGONS!\n"
               << "SURFACE " << surface_name_ << " IS SET INCORRECTLY\n";
        exit(1);
    }
    coefs_ = Surface::CalcSurfaceCoefficients(contour_);
    normal_ = Vector(coefs_.A_, coefs_.B_, coefs_.C_).Norm();
    x_bnd_ = Surface::GetBoundary(contour_, 'x');
    y_bnd_ = Surface::GetBoundary(contour_, 'y');
    z_bnd_ = Surface::GetBoundary(contour_, 'z');
}


Surface::Boundary Surface::GetBoundary(const std::vector<Vector>& ctr, char axis){
    std::vector<double> tmp;
    std::vector<double>::const_iterator res_min;
    std::vector<double>::const_iterator res_max;
    switch (axis) {
    case 'x':
        tmp.push_back(ctr[0].GetX());
        tmp.push_back(ctr[1].GetX());
        tmp.push_back(ctr[2].GetX());
        tmp.push_back(ctr[3].GetX());
        res_min = std::min_element(tmp.begin(), tmp.end());
        res_max = std::max_element(tmp.begin(), tmp.end());
        return {*res_min, *res_max};
    case 'y':
        tmp.push_back(ctr[0].GetY());
        tmp.push_back(ctr[1].GetY());
        tmp.push_back(ctr[2].GetY());
        tmp.push_back(ctr[3].GetY());
        res_min = std::min_element(tmp.begin(), tmp.end());
        res_max = std::max_element(tmp.begin(), tmp.end());
        return {*res_min, *res_max};
    case 'z':
        tmp.push_back(ctr[0].GetZ());
        tmp.push_back(ctr[1].GetZ());
        tmp.push_back(ctr[2].GetZ());
        tmp.push_back(ctr[3].GetZ());
        res_min = std::min_element(tmp.begin(), tmp.end());
        res_max = std::max_element(tmp.begin(), tmp.end());
        return {*res_min, *res_max};
    default:
        std::cerr << "Unknown axis!" << axis << "\n";
        exit(1);
    }
}

Surface::SurfaceCoeficients Surface::CalcSurfaceCoefficients(
                                            const std::vector<Vector> contour){
    double A = (contour[1].GetY() - contour[0].GetY())*(contour[2].GetZ()- contour[0].GetZ())
            - (contour[2].GetY() - contour[0].GetY())*(contour[1].GetZ()- contour[0].GetZ());
    double B =-(contour[1].GetX() - contour[0].GetX())*(contour[2].GetZ()- contour[0].GetZ())
            + (contour[2].GetX() - contour[0].GetX())*(contour[1].GetZ()- contour[0].GetZ());
    double C = (contour[1].GetX() - contour[0].GetX())*(contour[2].GetY() - contour[0].GetY())
            - (contour[2].GetX() - contour[0].GetX())*(contour[1].GetY() - contour[0].GetY());
    double D = -A*contour[0].GetX() - B*contour[0].GetY() - C*contour[0].GetZ();
    return {A, B, C, D};
}

const std::vector<Vector>& Surface::GetContour() const{return contour_;}
const Vector& Surface::GetNormal() const{return normal_;}
bool Surface::IsSaveStat() const{ return save_stat_;}
const Reflector* Surface::GetReflector() const {return reflector_.get();}
const std::string& Surface::GetName() const {return surface_name_;}
void Surface::SaveParticle(const Particle& pt){stat_.push_back(pt);}
const Surface::SurfaceCoeficients& Surface::GetSurfaceCoefficients() const {
    return coefs_;
}
const Surface::Boundary& Surface::GetXBnd() const {return x_bnd_;}
const Surface::Boundary& Surface::GetYBnd() const {return y_bnd_;}
const Surface::Boundary& Surface::GetZBnd() const {return z_bnd_;}



void Surface::SaveSurfaceParticles(std::ofstream& out) const{
    if(!stat_.empty()){
        out << "#POS_X\tPOS_Y\tPOS_Z\tVX\tVY\tVZ\tVolumeCount\tSurfaceCount\n";
        for(const auto& el : stat_){
            out << el.GetPosition() << "\t" << el.GetDirection()
                <<"\t" << el.GetVolCount() << "\t" << el.GetSurfCount() << "\n";
        }
    }
}

bool Surface::CheckIfPointOnSurface(const Vector& point) const{
    //TODO Make adequate check for polygon by counting rotation algo...
    if(x_bnd_.min_==x_bnd_.max_ &&
            point.GetY()<y_bnd_.max_ && point.GetY()>y_bnd_.min_ &&
            point.GetZ()<z_bnd_.max_ && point.GetZ()>z_bnd_.min_ ){
        return true;
    }
    if(y_bnd_.min_==y_bnd_.max_ &&
            point.GetX()<x_bnd_.max_ && point.GetX()>x_bnd_.min_ &&
            point.GetZ()<z_bnd_.max_ && point.GetZ()>z_bnd_.min_ ){
        return true;
    }
    if(z_bnd_.min_==z_bnd_.max_ &&
            point.GetX()<x_bnd_.max_ && point.GetX()>x_bnd_.min_ &&
            point.GetY()<y_bnd_.max_ && point.GetY()>y_bnd_.min_ ){
        return true;
    }
    return false;
}
