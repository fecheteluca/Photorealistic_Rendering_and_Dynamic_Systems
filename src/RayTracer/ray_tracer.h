#ifndef RAY_TRACER_H
#define RAY_TRACER_H

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0);
    double norm2() const;
    double norm() const;
    void normalize();

    double operator[](int i) const;
    double& operator[](int i);

private:
    double data[3];
};

// Non-member operator overloads for Vector
Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const double a, const Vector& b);
Vector operator*(const Vector& a, const double b);
Vector operator/(const Vector& a, const double b);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

// Stub declarations for additional classes
class Sphere {
public:
    // ... (Sphere implementation goes here)
};

class Ray {
public:
    // ... (Ray implementation goes here)
};

#endif // RAY_TRACER_H
