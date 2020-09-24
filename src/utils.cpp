#include <cmath>
#include "utils.hpp"


Vector::Vector(const Direction dir){
    switch (dir) {
    case Direction::POSITIVE_X:
        x_ = 1.0, y_ = 0.0, z_ = 0.0;
        break;
    case Direction::POSITIVE_Y:
        x_ = 0.0, y_ = 1.0, z_ = 0.0;
        break;
    case Direction::POSITIVE_Z:
        x_ = 0.0, y_ = 0.0, z_ = 1.0;
        break;
    case Direction::NEGATIVE_X:
        x_ = -1.0, y_ = 0.0, z_ = 0.0;
        break;
    case Direction::NEGATIVE_Y:
        x_ = 0.0, y_ = -1.0, z_ = 0.0;
        break;
    case Direction::NEGATIVE_Z:
        x_ = 0.0, y_ = 0.0, z_ = -1.0;
        break;
    default:
        stderr << "Setting for the given direcion is not added!\n";
        exit(1);
    }
}

double Vector::GetX() const {return x_;}
double Vector::GetY() const {return y_;}
double Vector::GetZ() const {return z_;}



Vector Vector::Dot(const Vector &rhs) const{
    return Vector(x_*rhs.GetX(), y_*rhs.GetY(), z_*rhs.GetZ());
}

Vector::Vector(const Vector &start_point, const Vector &end_point){
    return Vector(end_point.GetX() - start_point.GetX(),
                  end_point.GetY() - start_point.GetY(),
                  end_point.GetZ() - start_point.GetZ());
}

Vector Vector::Cross(const Vector& rhs) const{
    return Vector(y_*rhs.GetZ() - z_*rhs.GetY(),
                  z_*rhs.GetX() - x_*rhs.GetZ(),
                  x_*rhs.GetY() - y_*rhs.GetX());
}



double Vector::Length2(){return x_*x_ + y_*y_ + z_*z_;}
double Vector::Length(){return sqrt(Length2());}

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




