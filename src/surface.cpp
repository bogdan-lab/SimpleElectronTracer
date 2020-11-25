#include <utility>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>
#include <fmt/core.h>

#include "surface.hpp"

Surface::Surface(std::vector<Vec3>&& g_contour,
        std::unique_ptr<Reflector>&& g_reflector, std::ofstream&& out_file,
                  std::unique_ptr<char[]>&& buff):
    contour_(std::move(g_contour)),
    reflector_(std::move(g_reflector)),
    io_buffer_(std::move(buff)), output_file_(std::move(out_file))
{
    coefs_ = Surface::CalcSurfaceCoefficients(contour_);
    normal_ = Vec3(coefs_.A_, coefs_.B_, coefs_.C_).Norm();
    tri_areas_ = Surface::CalcTriangleAreas(contour_);
    mass_center_ = Surface::CalcCenterOfMass(contour_);
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
    return {x/contour.size(), y/contour.size(), z/contour.size()};
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
const Vec3& Surface::GetNormal() const{return normal_;}
bool Surface::IsSaveStat() const{ return output_file_.is_open();}
const Reflector* Surface::GetReflector() const {return reflector_.get();}
const Vec3& Surface::GetMassCenter() const{return mass_center_;}


const Surface::SurfaceCoeficients& Surface::GetSurfaceCoefficients() const {
    return coefs_;
}

void Surface::WriteFileHeader(){
    if(output_file_.is_open()){
        output_file_ << "#POS_X\tPOS_Y\tPOS_Z"
                     << "\tVX\tVY\tVZ\tVolumeCount\tSurfaceCount\n";
    }
}

void Surface::SaveParticle(Particle&& pt){
        output_file_ << fmt::format("{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:d}\t{:d}\n",
                                    pt.GetPosition().GetX(),
                                    pt.GetPosition().GetY(),
                                    pt.GetPosition().GetZ(),
                                    pt.GetDirection().GetX(),
                                    pt.GetDirection().GetY(),
                                    pt.GetDirection().GetZ(),
                                    pt.GetVolCount(),
                                    pt.GetSurfCount());

}

std::vector<Vec3> Surface::TranslateContourIntoBasis(
        const ONBasis_3x3 &basis) const{
    std::vector<Vec3> basis_contour;
    basis_contour.reserve(contour_.size()+1);
    for(size_t i=0; i<contour_.size(); i++){
        basis_contour.push_back(basis.FromOriginalCoorsToThis(contour_[i]));
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
    ONBasis_3x3 surf_basis(normal_);
    Vec3 basis_point = surf_basis.FromOriginalCoorsToThis(point);
    std::vector<Vec3> basis_contour = TranslateContourIntoBasis(surf_basis);
    if(!(basis_contour[0]==basis_contour[contour_.size()-1])){
        basis_contour.push_back(basis_contour[0]); //loop contour
    }
    //Counting rotations....
    std::vector<int> quarters = PrepareQuarterListForContour(basis_contour,
                                                             basis_point);
    int winding_num = CalculateWindingNumber(quarters, basis_contour,
                                             basis_point);
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
    auto defect = from_s_to_point.Dot(normal_);
    auto cos_alpha = normal_.Dot(pt_direction);
    while( defect<0){
        end = end - pt_direction.Times(defect/cos_alpha);
        from_s_to_point = Vec3(mass_center_, end);
        defect = from_s_to_point.Dot(normal_);
    }
}
