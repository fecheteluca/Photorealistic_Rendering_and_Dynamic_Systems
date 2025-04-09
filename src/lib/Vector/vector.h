#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <cmath>

class Vector {
    public:
        explicit Vector(
            double x = 0, 
            double y = 0, 
            double z = 0
        );

        void set_x(double x);
        void set_y(double y);
        void set_z(double z);

        double get_x();
        double get_y();
        double get_z();

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
    Vector operator*(const Vector &a, const Vector& b);
    Vector operator/(const Vector& a, const double b);
    double dot(const Vector& a, const Vector& b);
    Vector cross(const Vector& a, const Vector& b);

#endif // VECTOR_H