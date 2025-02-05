#include "vec_class.h"

vector3D::vector3D() : x(0), y(0), z(0) {}

vector3D::vector3D(double x_set, double y_set, double z_set) : x(x_set), y(y_set), z(z_set) {}

double vector3D::getX() const {
    return x;
}

double vector3D::getY() const {
    return y;
}

double vector3D::getZ() const {
    return z;
}

void vector3D::setX(double x) {
    this->x = x; // The "this->x" makes sure the x isn't an arbitrary parameter, but the particular x of the vector3D object
}

void vector3D::setY(double y) {
    this->y = y;
}

void vector3D::setZ(double z) {
    this->z = z;
}

// Could have avoided name shadowing by using different names, and thus avoided the need for `this->`
// void setZ(double newZ) { 
//     z = newZ; // No need for `this->` since there's no name conflict
// }

// this->z is equivalent to (*this).z

double vector3D::norm() const {
    return std::sqrt(this->dot(*this));
}

double vector3D::dot(const vector3D& other) const {
    return x * other.x + y * other.y + z * other.z;
}

vector3D vector3D::cross(const vector3D& other) const {
    return vector3D(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

void vector3D::normalize() {
    double mag = norm();
    if (mag != 0) {
        x /= mag;
        y /= mag;
        z /= mag;
    }
}


// Operator overloading:
vector3D vector3D::operator+(const vector3D& other) const {
    return vector3D(x + other.x, y + other.y, z + other.z);
}

vector3D vector3D::operator-(const vector3D& other) const {
    return vector3D(x - other.x, y - other.y, z - other.z);
}

vector3D vector3D::operator*(double scalar) const {
    return vector3D(scalar*x, scalar*y, scalar*z);
}

vector3D vector3D::operator/(double scalar) const {
    return vector3D(x/scalar, y/scalar, z/scalar);
}

// Left scalar multiplication
vector3D operator*(double scalar, const vector3D& vec) {
    return vec * scalar;
}