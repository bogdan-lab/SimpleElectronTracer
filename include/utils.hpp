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


class Vec3{
private:
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
public:
    Vec3(const double x, const double y, const double z):
        x_(x), y_(y), z_(z) {}
    Vec3(const Vec3& start_point, const Vec3& end_point);
    Vec3()=default;

    double GetX() const;
    double GetY() const;
    double GetZ() const;

    double Dot(const Vec3& rhs) const;
    Vec3 Cross(const Vec3& rhs) const;
    double Length2() const;
    double Length() const;
    Vec3 Norm() const;
    Vec3 Times(const double n) const;
};

Vec3 operator+(const Vec3& lhs, const Vec3& rhs);
Vec3 operator-(const Vec3& lhs, const Vec3& rhs);
bool operator==(const Vec3& lhs, const Vec3& rhs);
std::ostream& operator<<(std::ostream& out, const Vec3& vec);
double GetDistance(const Vec3& start, const Vec3& end);

class Basis_3x3{
private:
    std::vector<Vec3> m_;
public:
    Basis_3x3() = delete;
    Basis_3x3(const Vec3& i, const Vec3& j, const Vec3& k);
    Basis_3x3(const Vec3& given_z);
    Vec3 ApplyToVec(const Vec3& vec) const ;
    double GetDeterminant() const ;
    Basis_3x3 Transpose() const ;
    Basis_3x3 GetInverse() const;
    const std::vector<Vec3>& GetBasisCols() const;
};

Vec3 VerifyPointOnSurface(const Surface& s, const Vec3& point);
#endif //UTILS_HEADER
