#include <utility>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <fmt/core.h>
#include <omp.h>

#include "surface.hpp"

Surface::Surface(std::vector<Vec3>&& g_contour,
        std::unique_ptr<Reflector>&& g_reflector,
        bool g_save_stat_flag):
    contour_(std::move(g_contour)),
    reflector_(std::move(g_reflector)),
    save_stat_flag_(g_save_stat_flag)
{
    coefs_ = Surface::CalcSurfaceCoefficients(contour_);
    tri_areas_ = Surface::CalcTriangleAreas(contour_);
    mass_center_ = Surface::CalcCenterOfMass(contour_);
    surf_basis_ = ONBasis_3x3(Vec3(coefs_.A_, coefs_.B_, coefs_.C_).Norm());
    basis_contour_ = Surface::TranslateContourIntoBasis(surf_basis_, contour_);
}


std::vector<double> Surface::CalcTriangleAreas(const std::vector<Vec3> &contour){
    std::vector<double> areas;
    areas.reserve(contour.size()-2);
    for(size_t i=1; i<contour.size()-1; i++){
        areas.push_back(0.5*Vec3(contour[0], contour[i])
                .Cross(Vec3(contour[0], contour[i+1]))
                .Length());
    }
    return areas;
}

Vec3 Surface::CalcCenterOfMass(const std::vector<Vec3>& contour){
    double x = 0;
    double y = 0;
    double z = 0;
    for(size_t i=0; i<contour.size(); i++){
        x+=contour[i].GetX();
        y+=contour[i].GetY();
        z+=contour[i].GetZ();
    }
    return {x/static_cast<double>(contour.size()),
            y/static_cast<double>(contour.size()),
            z/static_cast<double>(contour.size())};
}

Vec3 Surface::GetRandomPointInContour(std::mt19937 &rng) const{
    std::discrete_distribution<size_t> dist(tri_areas_.begin(), tri_areas_.end());
    size_t tri_idx = dist(rng);
    std::uniform_real_distribution<double> dist_2(0.0, 1.0);
    double r1 = dist_2(rng);
    double r2 = dist_2(rng);
    return contour_[0].Times(1-sqrt(r1)) +
           contour_[tri_idx+1].Times(sqrt(r1)*(1-r2)) +
            contour_[tri_idx+2].Times(r2*sqrt(r1));
}

Surface::SurfaceCoeficients Surface::CalcSurfaceCoefficients(
                                            const std::vector<Vec3> contour){
    double A = (contour[1].GetY() - contour[0].GetY())*(contour[2].GetZ()- contour[0].GetZ())
            - (contour[2].GetY() - contour[0].GetY())*(contour[1].GetZ()- contour[0].GetZ());
    double B =-(contour[1].GetX() - contour[0].GetX())*(contour[2].GetZ()- contour[0].GetZ())
            + (contour[2].GetX() - contour[0].GetX())*(contour[1].GetZ()- contour[0].GetZ());
    double C = (contour[1].GetX() - contour[0].GetX())*(contour[2].GetY() - contour[0].GetY())
            - (contour[2].GetX() - contour[0].GetX())*(contour[1].GetY() - contour[0].GetY());
    double D = -A*contour[0].GetX() - B*contour[0].GetY() - C*contour[0].GetZ();
    return {A, B, C, D};
}

const std::vector<Vec3>& Surface::GetContour() const{return contour_;}
const Vec3& Surface::GetNormal() const{return surf_basis_.GetZVec();}
bool Surface::IsSaveStat() const{ return save_stat_flag_;}
const Reflector* Surface::GetReflector() const {return reflector_.get();}
const Vec3& Surface::GetMassCenter() const{return mass_center_;}


const Surface::SurfaceCoeficients& Surface::GetSurfaceCoefficients() const {
    return coefs_;
}

std::vector<Vec3> Surface::TranslateContourIntoBasis(
        const ONBasis_3x3 &basis, const std::vector<Vec3>& contour){
    std::vector<Vec3> basis_contour;
    basis_contour.reserve(contour.size()+1);
    for(size_t i=0; i<contour.size(); i++){
        basis_contour.push_back(basis.FromOriginalCoorsToThis(contour[i]));
    }
    return basis_contour;
}

bool Surface::CheckIfPointOnSurface(const Vec3& point) const{
    //Due to the finite double precision we still expect that given point
    //will be not directly on the surface
    //anyway we will transform its and surface coordinates into the basis
    //where Z is parallel to the normal and than compare X and Y coordinates
    //of the point and surface polygon in order to answer the question whether
    //point is on the surface
    Vec3 basis_point = surf_basis_.FromOriginalCoorsToThis(point);
    //Counting rotations....
    int winding_num = 0;
    for(size_t i=0; i<basis_contour_.size()-1; i++){
        winding_num += CalcWindChange(basis_contour_[i], basis_contour_[i+1],
                basis_point);
    }
    winding_num += CalcWindChange(basis_contour_.back(), basis_contour_.front(),
                                  basis_point);
    winding_num/=4;
    return !(winding_num%2==0);
}

std::optional<Vec3> Surface::GetCrossPoint(const Vec3& pos,
                                           const Vec3& dir) const {
    if((coefs_.A_*dir.GetX()
        + coefs_.B_*dir.GetY()
        + coefs_.C_*dir.GetZ()) == 0.0){
        //particle moves parallel to the surface
        return std::nullopt;
    }
    //Look at time needed to reach the surface
    double tmp_num = coefs_.A_*pos.GetX() + coefs_.B_*pos.GetY() +
                     coefs_.C_*pos.GetZ() + coefs_.D_;
    double tmp_den = coefs_.A_*dir.GetX()
            + coefs_.B_*dir.GetY() + coefs_.C_*dir.GetZ();
    double t = -1*tmp_num/tmp_den;
    if(t<=0){
        return std::nullopt;
    }
    //Here at least direction is correct --> check for boundaries
    Vec3 cross_point = {pos.GetX() + dir.GetX()*t,
                        pos.GetY() + dir.GetY()*t,
                        pos.GetZ() + dir.GetZ()*t};
    VerifyPointInVolume(pos, cross_point);
    if(CheckIfPointOnSurface(cross_point)){
        return cross_point;
    }
    return std::nullopt;
}

void Surface::VerifyPointInVolume(const Vec3& start, Vec3& end) const {
    /*!Function assumes that surface normal is directed inside the volume!*/
    Vec3 from_s_to_point(mass_center_, end);
    Vec3 pt_direction(start, end);
    pt_direction.Norm();
    auto defect = from_s_to_point.Dot(GetNormal());
    auto cos_alpha = GetNormal().Dot(pt_direction);
    while( defect<0){
        end = end - pt_direction.Times(defect/cos_alpha);
        from_s_to_point = Vec3(mass_center_, end);
        defect = from_s_to_point.Dot(GetNormal());
    }
}



int Surface::GetQuarter(const Vec3& point, const Vec3& node) const {
    if(node.GetX()>point.GetX() && node.GetY()>=point.GetY()){
        return 0;
    }
    if(node.GetX()<=point.GetX() && node.GetY()>point.GetY()){
        return 1;
    }
    if(node.GetX()<point.GetX() && node.GetY()<=point.GetY()){
        return 2;
    }
    if(node.GetX()>=point.GetX() && node.GetY()<point.GetY()){
        return 3;
    }
    fprintf(stderr, "INCORRECT QUARTER CALCULATION!");
    exit(1);
}


double Surface::GetOrientationWinding(const Vec3& point,
                         const Vec3& prev_node, const Vec3& next_node) const {
    double tmp_1 = (prev_node.GetX() - point.GetX())*
            (next_node.GetY() - point.GetY());
    double tmp_2 = (next_node.GetX() - point.GetX())*
            (prev_node.GetY() - point.GetY());
    return  tmp_1 - tmp_2;
}

int Surface::CalcWindChange(const Vec3& prev_node, const Vec3& next_node,
                            const Vec3& point) const {
    int quarter_prev = GetQuarter(point, prev_node);
    int quarter_next = GetQuarter(point, next_node);
    switch (quarter_next-quarter_prev){
    case 0:
        return 0;
    case 1:
    case -3:{
        return 1;
    }
    case -1:
    case 3:{
        return -1;
    };
    case 2:
    case -2:{
        double det = GetOrientationWinding(point, prev_node, next_node);
        return det>0 ? 2 : -2;
    }
    default:{
        fprintf(stderr, "Something went wrong in count algo");
        exit(1);
    }
    }
}
