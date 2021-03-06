﻿#include <cmath>
#include <vector>
#include <random>
#include <utility>

#include "math.hpp"
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

Vec3& Vec3::Norm() {
    double l = Length();
    x_ /= l;
    y_ /= l;
    z_ /= l;
    return *this;
}

Vec3 Vec3::Times(const double d) const{
    return Vec3(x_*d, y_*d, z_*d);
}

Vec3 Vec3::operator +(const Vec3& other) const {
    return Vec3(this->x_ + other.x_,
                this->y_ + other.y_,
                this->z_ + other.z_);
}

Vec3 Vec3::operator -(const Vec3& other) const {
    return Vec3(this->x_ - other.x_,
                this->y_ - other.y_,
                this->z_ - other.z_);
}

bool Vec3::operator ==(const Vec3& other) const {
    return this->x_==other.x_ && this->y_ == other.y_ && this->z_ == other.z_;
}

bool Vec3::operator !=(const Vec3& other) const{
    return this->x_!=other.x_ || this->y_ != other.y_ || this->z_ != other.z_;
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


ONBasis_3x3::ONBasis_3x3(Vec3 i, Vec3 j,Vec3 k):
    i_(std::move(i)), j_(std::move(j)), k_(std::move(k)){
    if(i.Dot(j)!=0 || j.Dot(k)!=0 || i.Dot(k)!=0){
        fprintf(stderr, "Basis is not orthogonal!\n");
        exit(1);
    }
    i_.Norm();
    j_.Norm();
    k_.Norm();
}


ONBasis_3x3::ONBasis_3x3(Vec3 new_z): k_(std::move(new_z)){
    Vec3 tmp_cross_x = k_.Cross(Vec3(1.0, 0.0, 0.0));
    Vec3 tmp_cross_y = k_.Cross(Vec3(0.0, 1.0, 0.0));
    j_ = tmp_cross_x.Length2()>tmp_cross_y.Length2() ?
                tmp_cross_x : tmp_cross_y;
    i_ = j_.Cross(k_);
    i_.Norm();
    j_.Norm();
    k_.Norm();
}

const Vec3& ONBasis_3x3::GetXVec() const{return i_;}
const Vec3& ONBasis_3x3::GetYVec() const{return j_;}
const Vec3& ONBasis_3x3::GetZVec() const{return k_;}


Vec3 ONBasis_3x3::ApplyToVec(const Vec3& vec) const {
    return i_.Times(vec.GetX()) +
           j_.Times(vec.GetY()) +
           k_.Times(vec.GetZ());
}

Vec3 ONBasis_3x3::FromOriginalCoorsToThis(const Vec3& vec) const {
    return {vec.Dot(i_), vec.Dot(j_), vec.Dot(k_)};
}

Vec3 ONBasis_3x3::FromThisCoorsToOriginal(const Vec3& vec) const{
    return i_.Times(vec.GetX()) +
           j_.Times(vec.GetY()) +
           k_.Times(vec.GetZ());
}



