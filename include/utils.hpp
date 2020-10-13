#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <vector>

class Surface;

struct Background{
    double sigma_; 	//crossection cm^2
    double T_; 	 	//temperature K
    double p_; 	 	//pressure Pa
};


class Vector{
private:
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
public:
    Vector(const double x, const double y, const double z):
        x_(x), y_(y), z_(z) {}
    Vector(const Vector& start_point, const Vector& end_point);
    Vector()=default;

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    double Dot(const Vector& rhs) const;
    Vector Cross(const Vector& rhs) const;
    double Length2() const;
    double Length() const;
    Vector Norm() const;
    Vector Times(const double n) const;
};

Vector operator+(const Vector& lhs, const Vector& rhs);
Vector operator-(const Vector& lhs, const Vector& rhs);
std::ostream& operator<<(std::ostream& out, const Vector& vec);
double GetDistance(const Vector& start, const Vector& end);
Vector VerifyPointOnSurface(const Surface& s, const Vector& point);
std::vector<Vector> GenerateONBasisByNewZ(const Vector& new_z);
Vector ApplyCoordinateTransition(const std::vector<Vector>& new_basis,
                                 const Vector& vec);
#endif //UTILS_HEADER
