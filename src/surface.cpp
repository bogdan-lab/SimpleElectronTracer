#include <utility>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>

#include "surface.hpp"

Surface::Surface(std::vector<Vec3> g_contour,
        std::unique_ptr<Reflector> g_reflector, FILE* out_file,
                 size_t dump_size):
    stat_(), contour_(std::move(g_contour)),
    reflector_(std::move(g_reflector)), output_file_(std::move(out_file))
{
    stat_.reserve(dump_size);
    coefs_ = Surface::CalcSurfaceCoefficients(contour_);
    normal_ = Vec3(coefs_.A_, coefs_.B_, coefs_.C_).Norm();
}

Surface::~Surface(){
    if(output_file_){
        SaveSurfaceParticles();
        fclose(output_file_);
    }
}

Vec3 Surface::GetPointOnSurface() const{
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
bool Surface::IsSaveStat() const{ return bool(output_file_);}
const Reflector* Surface::GetReflector() const {return reflector_.get();}

void Surface::SaveParticle(Particle&& pt){
    stat_.push_back(std::move(pt));
    if(stat_.size()==stat_.capacity()){
        SaveSurfaceParticles();
        stat_.clear();

    }
}

const Surface::SurfaceCoeficients& Surface::GetSurfaceCoefficients() const {
    return coefs_;
}

void Surface::WriteFileHeader() const{
    if(output_file_){
        fprintf(output_file_,
                "#POS_X\tPOS_Y\tPOS_Z\tVX\tVY\tVZ\tVolumeCount\tSurfaceCount\n");
    }
}

void Surface::SaveSurfaceParticles() const{
    for(const auto& pt : stat_){
            fprintf(output_file_,
                    "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%zu\t%zu\n",
                    pt.GetPosition().GetX(),
                    pt.GetPosition().GetY(),
                    pt.GetPosition().GetZ(),
                    pt.GetDirection().GetX(),
                    pt.GetDirection().GetY(),
                    pt.GetDirection().GetZ(),
                    pt.GetVolCount(),
                    pt.GetSurfCount());
        }
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
