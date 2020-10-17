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

Vec3 VerifyPointOnSurface(const Surface& s, const Vec3& point){
    //Stupid solution for now
    if(s.GetXBnd().max_==s.GetXBnd().min_){
        return Vec3(s.GetXBnd().min_, point.GetY(), point.GetZ());
    }
    else if (s.GetYBnd().min_ == s.GetYBnd().max_){
        return Vec3(point.GetX(), s.GetYBnd().max_, point.GetZ());
    }
    else if (s.GetZBnd().min_==s.GetZBnd().max_){
        return Vec3(point.GetX(), point.GetY(), s.GetZBnd().min_);
    }
    else{
        std::cerr << "Bended surfaces like " << s.GetName()
                  << "are not supported yet!\n";
        exit(1);
    }
}

Basis_3x3::Basis_3x3(Vec3 i, Vec3 j, Vec3 k){
    m_.push_back(std::move(i));
    m_.push_back(std::move(j));
    m_.push_back(std::move(k));
}


Basis_3x3::Basis_3x3(const Vec3& given_z){
    Vec3 tmp_cross_x = given_z.Cross(Vec3(1.0, 0.0, 0.0));
    Vec3 tmp_cross_y = given_z.Cross(Vec3(0.0, 1.0, 0.0));
    Vec3 new_y = tmp_cross_x.Length2()>tmp_cross_y.Length2() ?
                tmp_cross_x : tmp_cross_y;
    Vec3 new_x = new_y.Cross(given_z);
    m_.emplace_back(new_x);
    m_.emplace_back(new_y);
    m_.emplace_back(given_z);
}

const std::vector<Vec3>& Basis_3x3::GetBasisCols() const {return m_;}

Vec3 Basis_3x3::ApplyToVec(const Vec3& vec) const {
    return m_[0].Times(vec.GetX()) +
           m_[1].Times(vec.GetY()) +
           m_[2].Times(vec.GetZ());
}


double Basis_3x3::GetDeterminant() const{
    return m_[0].GetX()*(m_[1].GetY()*m_[2].GetZ() - m_[2].GetY()*m_[1].GetZ())
          -m_[1].GetX()*(m_[0].GetY()*m_[2].GetZ() - m_[2].GetY()*m_[0].GetZ())
          +m_[2].GetX()*(m_[0].GetY()*m_[1].GetZ() - m_[1].GetY()*m_[0].GetZ());
}

Basis_3x3 Basis_3x3::Transpose() const {
    return {Vec3(m_[0].GetX(), m_[1].GetX(), m_[2].GetX()),
            Vec3(m_[0].GetY(), m_[1].GetY(), m_[2].GetY()),
            Vec3(m_[0].GetZ(), m_[1].GetZ(), m_[2].GetZ())};
}

Basis_3x3 Basis_3x3::GetInverse() const {
    double det = GetDeterminant();
    if (det==0){
        fprintf(stderr, "Basis have zero deternminant!");
        exit(1);
    }
    double a11 = (m_[1].GetY()*m_[2].GetZ() - m_[2].GetY()*m_[1].GetZ())/det;
    double a12 = -(m_[0].GetY()*m_[2].GetZ() - m_[2].GetY()*m_[0].GetZ())/det;
    double a13 = (m_[0].GetY()*m_[1].GetZ() - m_[1].GetY()*m_[0].GetZ())/det;
    double a21 = -(m_[1].GetX()*m_[2].GetZ() - m_[2].GetX()*m_[1].GetZ())/det;
    double a22 = (m_[0].GetX()*m_[2].GetZ() - m_[2].GetX()*m_[0].GetZ())/det;
    double a23 = -(m_[0].GetX()*m_[1].GetZ() - m_[1].GetX()*m_[0].GetZ())/det;
    double a31 = (m_[1].GetX()*m_[2].GetY() - m_[2].GetX()*m_[1].GetY())/det;
    double a32 = -(m_[0].GetX()*m_[2].GetY() - m_[2].GetX()*m_[0].GetY())/det;
    double a33 = (m_[0].GetX()*m_[1].GetY() - m_[1].GetX()*m_[0].GetY())/det;
    Basis_3x3 tmp(Vec3(a11, a21, a31),
                   Vec3(a12, a22, a32),
                   Vec3(a13, a23, a33));
    return tmp.Transpose();
}


Basis_3x3 Basis_3x3::Norm() const {
    return {m_[0].Norm(), m_[1].Norm(), m_[2].Norm()};
}


