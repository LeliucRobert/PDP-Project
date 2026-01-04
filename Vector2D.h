#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <iosfwd>  // forward declare std::ostream
#include <cmath>

class Vector2D {
public:
    double x, y;

    Vector2D();
    Vector2D(double x_, double y_);

    Vector2D operator+(const Vector2D& other) const;
    Vector2D operator-(const Vector2D& other) const;
    Vector2D operator*(double scalar) const;

    Vector2D& operator+=(const Vector2D& other);
    Vector2D operator/(double scalar) const;


    double norm2() const;      // x^2 + y^2 (fast)
    double magnitude() const;  // sqrt(norm2())

    friend std::ostream& operator<<(std::ostream& os, const Vector2D& v);
};

#endif
