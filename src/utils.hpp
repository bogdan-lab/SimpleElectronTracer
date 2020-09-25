#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <iostream>
#include <vector>

struct Background{
    double sigma_; 	//crossection cm^2
    double T_; 	 	//temperature K
    double p_; 	 	//pressure Pa
};

//Just for now Don want to think about bended planes....
enum Direction{
    POSITIVE_X,
    POSITIVE_Y,
    POSITIVE_Z,
    NEGATIVE_X,
    NEGATIVE_Y,
    NEGATIVE_Z
};


class Vector{
private:
    double x_;
    double y_;
    double z_;
public:
    Vector(const double x, const double y, const double z):
        x_(x), y_(y), z_(z) {}
    Vector(const Direction dir);
    Vector(const Vector& start_point, const Vector& end_point);

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

#endif //UTILS_HEADER
