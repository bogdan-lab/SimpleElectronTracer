
#include <utility>
#include <vector>
#include <iostream>
#include <list>

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


Surface::Boundary GetBoundary(const std::vector<Vector>& ctr, char axis){
    switch (axis) {
    case 'x':
        std::vector<double> tmp = {ctr[0].GetX(), ctr[1].GetX(),
                          ctr[2].GetX(), ctr[3].GetX()};
        return {std::min(tmp.begin(), tmp.end()),
                std::max(tmp.begin(), tmp.end())};
        break;
    case 'y':
        std::vector<double> tmp = {ctr[0].GetY(), ctr[1].GetY(),
                          ctr[2].GetY(), ctr[3].GetY()};
        return {std::min(tmp.begin(), tmp.end()),
                std::max(tmp.begin(), tmp.end())};
        break;
    case 'z':
        std::vector<double> tmp = {ctr[0].GetZ(), ctr[1].GetZ(),
                          ctr[2].GetZ(), ctr[3].GetZ()};
        return {std::min(tmp.begin(), tmp.end()),
                std::max(tmp.begin(), tmp.end())};
        break;
    default:
        std::cerr << "Unknown axis!" << axis << "\n";
        exit(1);
    }
}

Surface::SurfaceCoeficients Surface::CalcSurfaceCoefficients(
                                            const std::vector<Point> contour){
    double A = (contour[1].y_ - contour_[0].y_)*(contour[2].z_ - contour_[0].z_)
            - (contour[2].y_ - contour_[0].y_)*(contour[1].z_ - contour_[0].z_);
    double B =-(contour[1].x_ - contour_[0].x_)*(contour[2].z_ - contour_[0].z_)
            + (contour[2].x_ - contour_[0].x_)*(contour[1].z_ - contour_[0].z_);
    double C = (contour[1].x_ - contour_[0].x_)*(contour[2].y_ - contour_[0].y_)
            - (contour[2].x_ - contour_[0].x_)*(contour[1].y_ - contour_[0].y_);
    double D = -A*contour[0].x_ - B*contour[0].y_ - C*contour[0].z_;
    return {A, B, C, D};
}

const std::vector<Vector>& Surface::GetContour() const{return contour_;}
const Vector& Surface::GetNormal() const{return normal_;}
bool Surface::IsSaveStat() const{ return save_stat_;}
const Reflector* Surface::GetReflector() const {return reflector_.get();}
const std::string Surface::GetName() const {return surface_name_;}
void Surface::SaveParticle(const Particle& pt){stat_.push_back(pt);}


void Surface::SaveSurfaceParticles(std::ofstream& out) const{
    for(const auto el : stat_){
        out << el.GetPosition() << "\t" << el.GetDirection()
            <<"\t" << el.GetVolCount() << "\t" << el.GetSurfCount() << "\n";
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
