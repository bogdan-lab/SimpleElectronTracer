
#include <utility>
#include <vector>
#include <iostream>
#include <list>
#include <algorithm>

#include "surface.hpp"

Surface::Surface(std::string g_name, std::vector<Vec3> g_contour,
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
    normal_ = Vec3(coefs_.A_, coefs_.B_, coefs_.C_).Norm();
}

Vec3 Surface::GetPointOnSurface() const{
    if(coefs_.C_!=0){
        return {0.0, 0.0, -coefs_.D_/coefs_.C_};
    } else if (coefs_.B_!=0){
        return {0.0, -coefs_.D_/coefs_.B_, 0.0};
    } else if (coefs_.A_!=0){
        return {-coefs_.D_/coefs_.A_, 0.0, 0.0};
    } else {
        fprintf(stderr, "Surface %s has all coefficients equal to zero!!",
                surface_name_.c_str());
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
bool Surface::IsSaveStat() const{ return save_stat_;}
const Reflector* Surface::GetReflector() const {return reflector_.get();}
const std::string& Surface::GetName() const {return surface_name_;}
void Surface::SaveParticle(const Particle& pt){stat_.push_back(pt);}
const Surface::SurfaceCoeficients& Surface::GetSurfaceCoefficients() const {
    return coefs_;
}


void Surface::SaveSurfaceParticles() const{
    if(!stat_.empty()){
        using FileHolder = std::unique_ptr<FILE, int(*)(FILE*)>;
        FileHolder output_file(fopen(surface_name_.c_str(), "w"), fclose);
        if(!output_file){
            fprintf(stderr, "Could not create file %s\n", surface_name_.c_str());
            exit(1);
        }
        fprintf(output_file.get(),
                "#POS_X\tPOS_Y\tPOS_Z\tVX\tVY\tVZ\tVolumeCount\tSurfaceCount\n");
        for(const auto& pt : stat_){
            fprintf(output_file.get(),
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
}

bool Surface::CheckIfPointOnSurface(const Vec3& point) const{
    //Due to the finite double precision we still expect that given point
    //will be not directly on the surface
    //anyway we will transform its and surface coordinates into the basis
    //where Z is parallel to the normal and than compare X and Y coordinates
    //of the point and surface polygon in order to answer the question whether
    //point is on the surface
    Basis_3x3 surf_basis(normal_);
    surf_basis = surf_basis.Norm();
    Vec3 basis_point = surf_basis.GetCoordinatesInThis(point);
    std::vector<Vec3> basis_contour;
    basis_contour.reserve(contour_.size()+1);
    for(size_t i=0; i<contour_.size(); i++){
        basis_contour.push_back(surf_basis.GetCoordinatesInThis(contour_[i]));
    }
    if(basis_contour[0]!=basis_contour[contour_.size()-1]){
        basis_contour.push_back(basis_contour[0]); //loop contour
    }
    //Counting rotations....
    //prepare quater information for each point
    std::vector<int> quarters(basis_contour.size(), 0);
    for(size_t i=0; i<basis_contour.size(); i++){
        if(basis_contour[i].GetX()>basis_point.GetX() &&
                basis_contour[i].GetY()>=basis_point.GetY()){
            quarters[i] = 0;
        }
        if(basis_contour[i].GetX()<=basis_point.GetX() &&
                basis_contour[i].GetY()>basis_point.GetY()){
            quarters[i] = 1;
        }
        if(basis_contour[i].GetX()<basis_point.GetX() &&
                basis_contour[i].GetY()<=basis_point.GetY()){
            quarters[i] = 2;
        }
        if(basis_contour[i].GetX()>=basis_point.GetX() &&
                basis_contour[i].GetY()<basis_point.GetY()){
            quarters[i] = 3;
        }
        fprintf(stderr, "INCORRECT QUARTER CALCULATION!");
        exit(1);
    }
    auto get_orient = [basis_point](const Vec3& nod_prev, const Vec3& nod_next){
        double tmp_1 = (nod_prev.GetX() - basis_point.GetX())*
                       (nod_next.GetY() - basis_point.GetY());
        double tmp_2 = (nod_next.GetX() - basis_point.GetX())*
                       (nod_prev.GetY() - basis_point.GetY());
        return  tmp_1 - tmp_2;
    };
    //Counting winding number
    int winding_num = 0;
    for(size_t i=0; i<quarters.size()-1; i++){
        switch (q[i+1]-q[i]){
            case 1:
            case -3:
                winding_num++;
                break;
            case -1:
            case 3:
                winding_num--;
                break;
            case 2:
            case -2:
                double det = get_orient(basis_contour[i], basis_contour[i+1]);
                winding_num += 2*abs(det)/det;
                break;
            default:
                fprintf(stderr, "Something went wrong in count algo");
                exit(1);
        }
    }
    winding_num/=4;
    return winding_num%2!=0;
}
