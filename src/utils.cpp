#include <cmath>
#include <vector>
#include <random>

#include "utils.hpp"
#include "surface.hpp"

double Vector::GetX() const {return x_;}
double Vector::GetY() const {return y_;}
double Vector::GetZ() const {return z_;}



double Vector::Dot(const Vector &rhs) const{
    return x_*rhs.GetX() + y_*rhs.GetY() + z_*rhs.GetZ();
}

Vector::Vector(const Vector &start_point, const Vector &end_point){
    x_ = end_point.GetX() - start_point.GetX();
    y_ = end_point.GetY() - start_point.GetY();
    z_ = end_point.GetZ() - start_point.GetZ();
}

Vector Vector::Cross(const Vector& rhs) const{
    return Vector(y_*rhs.GetZ() - z_*rhs.GetY(),
                  z_*rhs.GetX() - x_*rhs.GetZ(),
                  x_*rhs.GetY() - y_*rhs.GetX());
}



double Vector::Length2() const {return x_*x_ + y_*y_ + z_*z_;}
double Vector::Length() const {return sqrt(Length2());}

Vector Vector::Norm() const {
    double l = Length();
    return Vector(x_/l, y_/l, z_/l);
}

Vector Vector::Times(const double n) const{return Vector(n*x_, n*y_, n*z_);}

Vector operator+(const Vector& lhs, const Vector& rhs){
    return Vector(lhs.GetX() + rhs.GetX(),
                  lhs.GetY() + rhs.GetY(),
                  lhs.GetZ() + rhs.GetZ());
}

Vector operator-(const Vector& lhs, const Vector& rhs){
    return Vector(lhs.GetX() - rhs.GetX(),
                  lhs.GetY() - rhs.GetY(),
                  lhs.GetZ() - rhs.GetZ());
}


std::ostream& operator<<(std::ostream& out, const Vector& vec){
    out << vec.GetX() << "\t" << vec.GetY() << "\t" << vec.GetZ();
    return out;
}

double GetDistance(const Vector& start, const Vector& end){
    return Vector(start, end).Length();
}

Vector VerifyPointOnSurface(const Surface& s, const Vector& point){
    //Stupid solution for now
    if(s.GetXBnd().max_==s.GetXBnd().min_){
        return Vector(s.GetXBnd().min_, point.GetY(), point.GetZ());
    }
    else if (s.GetYBnd().min_ == s.GetYBnd().max_){
        return Vector(point.GetX(), s.GetYBnd().max_, point.GetZ());
    }
    else if (s.GetZBnd().min_==s.GetZBnd().max_){
        return Vector(point.GetX(), point.GetY(), s.GetZBnd().min_);
    }
    else{
        std::cerr << "Bended surfaces like " << s.GetName()
                  << "are not supported yet!\n";
        exit(1);
    }
}

std::vector<Vector> GenerateONBasisByNewZ(const Vector &new_z){
    std::vector<Vector> basis(3);
    basis[2] = new_z.Norm();
    Vector tmp_cross_x = new_z.Cross(Vector(1.0, 0.0, 0.0));
    Vector tmp_cross_y = new_z.Cross(Vector(0.0, 1.0, 0.0));
    Vector new_y = tmp_cross_x.Length2()>tmp_cross_y.Length2() ?
                tmp_cross_x : tmp_cross_y;
    basis[1] = new_y.Norm();
    basis[0] = basis[1].Cross(basis[2]);
    return basis;
}


Vector ApplyCoordinateTransition(const std::vector<Vector>& new_basis,
                                 const Vector& vec){
    return new_basis[0].Times(vec.GetX()) +
           new_basis[1].Times(vec.GetY()) +
           new_basis[2].Times(vec.GetZ());
}
