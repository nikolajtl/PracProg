#ifndef VEC_H
#define VEC_H

#include <cmath>  // For mathematical functions
#include <numbers> // For constants like std::numbers::pi

class vector3D{    
    public:
        double x, y, z; // Needed to allow deafult values

        vector3D();
        vector3D(double x_set, double y_set, double z_set); //NB: 'Should' be called i_set, but can demonstrate that the names need not match

        double getX() const;
        double getY() const;
        double getZ() const;

        void setX(double x);
        void setY(double y);
        void setZ(double z);

        double norm() const;

        double dot(const vector3D& other) const;

        vector3D cross(const vector3D& other) const;

        void normalize();

        // Operator overloading:
        vector3D operator+(const vector3D& other) const;
        vector3D operator-(const vector3D& other) const;
        vector3D operator*(double scalar) const;
        vector3D operator/(double scalar) const;

};

vector3D operator*(double scalar, const vector3D& vec);

#endif // VEC_H