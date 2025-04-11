#include "vector.h"

// Vector class definitions
Vector::Vector(
    double x, 
    double y, 
    double z
) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

void Vector::set_x(double x) {
    data[0] = x;
}

void Vector::set_y(double y) {
    data[1] = y;
}

void Vector::set_z(double z) {
    data[2] = z;
}

double Vector::get_x() {
    return data[0];
}

double Vector::get_y() {
    return data[1];
}

double Vector::get_z() {
    return data[2];
}

double Vector::norm2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
}

double Vector::norm() const {
    return sqrt(norm2());
}

void Vector::normalize() {
    double n = norm();
    data[0] /= n;
    data[1] /= n;
    data[2] /= n;
}

double Vector::operator[](int i) const {
    return data[i];
}

double& Vector::operator[](int i) {
    return data[i];
}

bool Vector::operator==(const Vector& a) {
    return (data[0] == a[0]) && (data[1] == a[1]) && (data[2] == a[2]);
}

// Non-member operator overloads for Vector
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector& a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(const Vector &a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1],
                  a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
}
