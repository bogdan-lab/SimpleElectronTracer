#ifndef UTILS_HEADER
#define UTILS_HEADER

#include <iostream>
#include <vector>

struct Point{
    double x_;
    double y_;
    double z_;
};

struct Velocity{
    double vx_;
    double vy_;
    double vz_;
};

enum Direction{
    POSITIVE_X,
    POSITIVE_Y,
    POSITIVE_Z,
    NEGATIVE_X,
    NEGATIVE_Y,
    NEGATIVE_Z
};

std::vector<double> GetDirectionVector(const Direction dir){
    switch (dir) {
    case Direction::POSITIVE_X:
        return {1.0, 0.0, 0.0};
    case Direction::POSITIVE_Y:
        return {0.0, 1.0, 0.0};
    case Direction::POSITIVE_Z:
        return {0.0, 0.0, 1.0};
    case Direction::NEGATIVE_X:
        return {-1.0, 0.0, 0.0};
    case Direction::NEGATIVE_Y:
        return {0.0, -1.0, 0.0};
    case Direction::NEGATIVE_Z:
        return {0.0, 0.0, -1.0};
    }
    stderr << "Setting for the given direcion is not added!\n";
    exit(1);
}

#endif //UTILS_HEADER
