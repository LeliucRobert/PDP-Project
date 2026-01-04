#pragma once
#include <array>
#include <ostream>
#include "Vector2D.h"

class Body {
public:
    double mass;        // use double for physics (gravitational effects)
    Vector2D position;
    Vector2D velocity;

    Body() : mass(1.0), position(), velocity() {}

    Body(double mass_, Vector2D position_, Vector2D velocity_)
        : mass(mass_), position(position_), velocity(velocity_) {}

    std::array<double, 5> serialize() const;
    static Body deserialize(const std::array<double, 5>& buffer);

    friend std::ostream& operator<<(std::ostream& os, const Body& b) {
        os << "mass: " << b.mass
           << " position: (" << b.position.x << ", " << b.position.y << ")"
           << " velocity: (" << b.velocity.x << ", " << b.velocity.y << ")";
        return os;
    }
};
