#pragma once
#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <algorithm>

#include "vector.h"

struct PolygonMoments {
    double area;
    double signed_area;
    Vector centroid;
};

class Polygon {  
public:
    std::vector<Vector> vertices;

    explicit Polygon(std::vector<Vector> aux_vertices = std::vector<Vector>());  
};

PolygonMoments compute_moments(const Polygon& polygon);  

Polygon polygon_clipping(const Polygon &subject_polygon, const Polygon &clip_polygon);

std::vector<Polygon> voronoi_ple(const std::vector<Vector> &sites);

std::vector<Polygon> weighted_voronoi_ple(const std::vector<Vector> &sites, const std::vector<double> &weights);

// Helper: Clip a polygon by a disk centered at disk_center of radius disk_radius
Polygon clip_by_disk(const Polygon& polygon, const Vector& disk_center, double disk_radius);

#endif // POLYGON_H