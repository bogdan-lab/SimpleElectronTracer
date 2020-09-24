
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
    if(contour_.size()<3){
        stderr << "Not enough points for creating surface "
               << surface_name_ << "\n";
        exit(1);
    }
    coefs_ = Surface::CalcSurfaceCoefficients(contour_);
    //TODO add calculation of the normal to the surface
}

const Surface::SurfaceCoeficients& Surface::CalcSurfaceCoefficients(
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
bool Surface::IsSaveStat() const{ return save_stat_;}


void Surface::SaveParticle(const Particle& pt){stat_.push_back(pt);}


void Surface::SaveSurfaceParticles(std::ofstream& out) const{
    for(const auto el : stat_){
        out << el.GetPosition() << "\t" << el.GetDirection()
            <<"\t" << el.GetVolCount() << "\t" << el.GetSurfCount() << "\n";
    }
}
