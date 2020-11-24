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
}

Vec3 Surface::GetPointOnSurface() const{
    //TODO if you save surface bassi you can get point on surface by it
    if(coefs_.C_!=0){
        return {0.0, 0.0, -coefs_.D_/coefs_.C_};
    } else if (coefs_.B_!=0){
        return {0.0, -coefs_.D_/coefs_.B_, 0.0};
    } else if (coefs_.A_!=0){
        return {-coefs_.D_/coefs_.A_, 0.0, 0.0};
    } else {
        fprintf(stderr, "Surface has all coefficients equal to zero!!");
        exit(1);
    }
}

Vec3 Surface::GetRandomPointInContour(std::mt19937 &rng) const{
    //TODO finish this function
    //Plan --> roll two numbers in range of bounding rectangle of the contour in surface basis
    return {};
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
    Vec3 p_on_s = GetPointOnSurface();
    Vec3 from_s_to_point(p_on_s, end);
    Vec3 pt_direction(start, end);
    pt_direction.Norm();
    auto defect = from_s_to_point.Dot(normal_);
    auto cos_alpha = normal_.Dot(pt_direction);
    while( defect<0){
        end = end - pt_direction.Times(defect/cos_alpha);
        from_s_to_point = Vec3(p_on_s, end);
        defect = from_s_to_point.Dot(normal_);
    }
}
