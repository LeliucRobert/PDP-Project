#include "Vector2D.h"
#include <ostream>

Vector2D::Vector2D() : x(0.0), y(0.0) {}
Vector2D::Vector2D(double x_, double y_) : x(x_), y(y_) {}

Vector2D Vector2D::operator+(const Vector2D& other) const {
    return Vector2D(x + other.x, y + other.y);
}

Vector2D Vector2D::operator-(const Vector2D& other) const {
    return Vector2D(x - other.x, y - other.y);
}

Vector2D Vector2D::operator*(double scalar) const {
    return Vector2D(x * scalar, y * scalar);
}

Vector2D& Vector2D::operator+=(const Vector2D& other) {
    x += other.x;
    y += other.y;
    return *this;
}

double Vector2D::norm2() const {
    return x * x + y * y;
}

double Vector2D::magnitude() const {
    return std::sqrt(norm2());
}

std::ostream& operator<<(std::ostream& os, const Vector2D& v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}
