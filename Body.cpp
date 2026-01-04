#include "Body.h"

std::array<double, 5> Body::serialize() const {
    return { mass, position.x, position.y, velocity.x, velocity.y };
}

Body Body::deserialize(const std::array<double, 5>& buffer) {
    Body b;
    b.mass = buffer[0];
    b.position.x = buffer[1];
    b.position.y = buffer[2];
    b.velocity.x = buffer[3];
    b.velocity.y = buffer[4];
    return b;
}
