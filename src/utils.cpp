#include <cmath>
#include <vector>
#include <random>
#include <utility>

#include "utils.hpp"
#include "surface.hpp"

double Vec3::GetX() const {return x_;}
double Vec3::GetY() const {return y_;}
double Vec3::GetZ() const {return z_;}



double Vec3::Dot(const Vec3 &rhs) const{
    return x_*rhs.GetX() + y_*rhs.GetY() + z_*rhs.GetZ();
}

Vec3::Vec3(const Vec3 &start_point, const Vec3 &end_point){
    x_ = end_point.GetX() - start_point.GetX();
    y_ = end_point.GetY() - start_point.GetY();
    z_ = end_point.GetZ() - start_point.GetZ();
}

Vec3 Vec3::Cross(const Vec3& rhs) const{
    return Vec3(y_*rhs.GetZ() - z_*rhs.GetY(),
                  z_*rhs.GetX() - x_*rhs.GetZ(),
                  x_*rhs.GetY() - y_*rhs.GetX());
}



double Vec3::Length2() const {return x_*x_ + y_*y_ + z_*z_;}
double Vec3::Length() const {return sqrt(Length2());}

Vec3 Vec3::Norm() const {
    double l = Length();
    return Vec3(x_/l, y_/l, z_/l);
}

Vec3 Vec3::Times(const double n) const{return Vec3(n*x_, n*y_, n*z_);}

Vec3 operator+(const Vec3& lhs, const Vec3& rhs){
    return Vec3(lhs.GetX() + rhs.GetX(),
                  lhs.GetY() + rhs.GetY(),
                  lhs.GetZ() + rhs.GetZ());
}

Vec3 operator-(const Vec3& lhs, const Vec3& rhs){
    return Vec3(lhs.GetX() - rhs.GetX(),
                  lhs.GetY() - rhs.GetY(),
                  lhs.GetZ() - rhs.GetZ());
}

bool operator==(const Vec3& lhs, const Vec3& rhs){
    return lhs.GetX()==rhs.GetX() &&
           lhs.GetY()==rhs.GetY() &&
           lhs.GetZ()==rhs.GetZ();
}

std::ostream& operator<<(std::ostream& out, const Vec3& vec){
    out << vec.GetX() << "\t" << vec.GetY() << "\t" << vec.GetZ();
    return out;
}

double Vec3::GetDistance(const Vec3& end) const {
    return sqrt((x_-end.GetX())*(x_-end.GetX()) +
                (y_-end.GetY())*(y_-end.GetY()) +
                (z_-end.GetZ())*(z_-end.GetZ()));
}

void VerifyPointInVolume(const Surface& s, Vec3& point){
    /*!Function assumes that surface normal is directed inside the volume!*/
    Vec3 p_on_s = s.GetPointOnSurface();
    Vec3 from_s_to_point(p_on_s, point);
    while(from_s_to_point.Dot(s.GetNormal())<0){
        point = point - s.GetNormal().Times(from_s_to_point.Dot(s.GetNormal()));
        from_s_to_point = Vec3(p_on_s, point);
    }
}

ONBasis_3x3::ONBasis_3x3(const Vec3& i, const Vec3& j, const Vec3& k){
    if(i.Dot(j)==0 && j.Dot(k)==0 && i.Dot(k)==0){
        m_.push_back(i.Norm());
        m_.push_back(j.Norm());
        m_.push_back(k.Norm());
    }
    else{
        fprintf(stderr, "Basis is not orthogonal!\n");
        exit(1);
    }
}


ONBasis_3x3::ONBasis_3x3(const Vec3& given_z){
    m_.reserve(3);
    Vec3 new_z = given_z.Norm();
    Vec3 tmp_cross_x = new_z.Cross(Vec3(1.0, 0.0, 0.0));
    Vec3 tmp_cross_y = new_z.Cross(Vec3(0.0, 1.0, 0.0));
    Vec3 new_y = tmp_cross_x.Length2()>tmp_cross_y.Length2() ?
                tmp_cross_x : tmp_cross_y;
    Vec3 new_x = new_y.Cross(new_z);
    m_.push_back(new_x.Norm());
    m_.push_back(new_y.Norm());
    m_.push_back(std::move(new_z));
}

const Vec3& ONBasis_3x3::GetXVec() const{return m_[0];}
const Vec3& ONBasis_3x3::GetYVec() const{return m_[1];}
const Vec3& ONBasis_3x3::GetZVec() const{return m_[2];}


Vec3 ONBasis_3x3::ApplyToVec(const Vec3& vec) const {
    return m_[0].Times(vec.GetX()) +
           m_[1].Times(vec.GetY()) +
           m_[2].Times(vec.GetZ());
}

Vec3 ONBasis_3x3::FromOriginalCoorsToThis(const Vec3& vec) const {
    return {vec.Dot(m_[0]), vec.Dot(m_[1]), vec.Dot(m_[2])};
}

Vec3 ONBasis_3x3::FromThisCoorsToOriginal(const Vec3& vec) const{
    return m_[0].Times(vec.GetX()) +
           m_[1].Times(vec.GetY()) +
           m_[2].Times(vec.GetZ());
}

int GetQuarter(const Vec3& point, const Vec3& node){
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


std::vector<int> PrepareQuarterListForContour(const std::vector<Vec3> &contour,
                                              const Vec3 &point){
    std::vector<int> quarters(contour.size(), 0);
    for(size_t i=0; i<contour.size(); i++){
        quarters[i] = GetQuarter(point, contour[i]);
    }
    return quarters;
}

int CalculateWindingNumber(const std::vector<int>& quarters,
                           const std::vector<Vec3> contour,
                           const Vec3& point){
    auto get_orient = [point](const Vec3& nod_prev, const Vec3& nod_next){
        double tmp_1 = (nod_prev.GetX() - point.GetX())*
                       (nod_next.GetY() - point.GetY());
        double tmp_2 = (nod_next.GetX() - point.GetX())*
                       (nod_prev.GetY() - point.GetY());
        return  tmp_1 - tmp_2;
    };
    //Counting winding number
    int winding_num = 0;
    for(size_t i=0; i<quarters.size()-1; i++){
        switch (quarters[i+1]-quarters[i]){
            case 0:
                break;
            case 1:
            case -3:{
                winding_num++;
            }break;
            case -1:
            case 3:{
                winding_num--;
            }break;
            case 2:
            case -2:{
                double det = get_orient(contour[i], contour[i+1]);
                winding_num += 2*abs(det)/det;
            }break;
            default:{
                fprintf(stderr, "Something went wrong in count algo");
                exit(1);
            }
        }
    }
    winding_num/=4;
    return winding_num;
}
