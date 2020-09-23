#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <iostream>
#include <vector>

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

    Vector Dot(const Vector& rhs) const;
    Vector Cross(const Vector& rhs) const;

};

Vector operator+(const Vector& lhs, const Vector& rhs);
Vector operator-(const Vector& lhs, const Vector& rhs);


#endif //UTILS_HEADER
