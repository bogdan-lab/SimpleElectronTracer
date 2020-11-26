#ifndef MATH_HPP
#define MATH_HPP

#include <iostream>
#include <vector>
#include <memory>
#include <math.h>

class Surface;

struct Background{
    double sigma_; 	//crossection cm^2
    double T_; 	 	//temperature K
    double p_; 	 	//pressure Pa
};


class Vec3{
private:
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
public:
    Vec3(const double x, const double y, const double z):
        x_(x), y_(y), z_(z) {}
    Vec3(const std::vector<double>& vec):
        x_(vec[0]), y_(vec[1]), z_(vec[2]) {}
    Vec3(const Vec3& start_point, const Vec3& end_point);
    Vec3()=default;

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    double Dot(const Vec3& rhs) const;
    Vec3 Cross(const Vec3& rhs) const;
    double Length2() const;
    double Length() const;
    Vec3& Norm();
    Vec3 Times(const double d) const ;
    double GetDistance(const Vec3& end) const;

    Vec3 operator +(const Vec3& other) const;
    Vec3 operator -(const Vec3& other) const;
    bool operator ==(const Vec3& other) const;
    bool operator !=(const Vec3& other) const;
};
std::ostream& operator <<(std::ostream& out, const Vec3& vec);

class ONBasis_3x3{
private:
    Vec3 i_;
    Vec3 j_;
    Vec3 k_;
public:
    ONBasis_3x3(): i_(Vec3(1.0, 0.0, 0.0)),
                   j_(Vec3(0.0, 1.0, 0.0)),
                   k_(Vec3(0.0, 0.0, 1.0)) {}
    ONBasis_3x3(Vec3 i, Vec3 j, Vec3 k);
    explicit ONBasis_3x3(Vec3 given_z);
    Vec3 ApplyToVec(const Vec3& vec) const ;
    Vec3 FromOriginalCoorsToThis(const Vec3& vec) const;
    Vec3 FromThisCoorsToOriginal(const Vec3& vec) const;
    const Vec3& GetXVec() const;
    const Vec3& GetYVec() const;
    const Vec3& GetZVec() const;
};

#endif //MATH_HEADER
