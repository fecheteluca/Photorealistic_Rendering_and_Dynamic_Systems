#include "polygon.h"

// Polygon class method defintions
Polygon::Polygon(std::vector<Vector> aux_vertices) {
    vertices = aux_vertices;
}

// Helper for computing different values from the vertices of a polygon
PolygonMoments compute_moments(const Polygon& polygon) {
    double area = 0.0;
    double centroid_x = 0.0;
    double centroid_y = 0.0;
    
    for (int i = 0; i < polygon.vertices.size(); i++) {
        const Vector& current_point = polygon.vertices[i];
        const Vector& next_point    = polygon.vertices[(i + 1) % polygon.vertices.size()];
        double d = current_point[0] * next_point[1] - next_point[0] * current_point[1];  
        
        area += d;
        centroid_x += (current_point[0] + next_point[0]) * d;
        centroid_y += (current_point[1] + next_point[1]) * d;
    }
    
    double signed_area = area * 0.5;
    area = std::abs(area) * 0.5;
    
    return {area, signed_area, Vector(centroid_x / (6 * area), centroid_y/(6 * area))};
}

// Polygon Clipping (Sutherland–Hodgman)
Vector intersect(const Vector &point_A, const Vector &point_B, Vector &point_u, Vector &point_v) {
    Vector normal(point_v.get_y() - point_u.get_y(), point_u.get_x() - point_v.get_x()); 
    double distance = dot(point_u - point_A, normal) / dot(point_B - point_A, normal);
    if (!(distance >= 0 && distance <= 1)) return point_A;
    return point_A + (point_B - point_A) * distance;
}

bool inside(const Vector &point_P, Vector &point_u, Vector &point_v) {
    Vector normal(point_v.get_y() - point_u.get_y(), point_u.get_x() - point_v.get_x());
    return dot(point_P - point_u, normal) <= 0;
}

Polygon polygon_clipping(const Polygon &subject_polygon, const Polygon &clip_polygon) {
    Polygon clipped_polygon = subject_polygon;

    for (int i = 0; i < clip_polygon.vertices.size(); i++) {
        Polygon out_polygon;

        // Due to the fact that we need to iterate over the edges of the clip polygon
        // (based on the pseudo-code at page 88), an edge being represented by 2 vertices,
        // we extract consecutive vertices of the clip polygon.
        Vector point_u = clip_polygon.vertices[i];
        Vector point_v = clip_polygon.vertices[(i + 1) % clip_polygon.vertices.size()];

        for (int j = 0; j < clipped_polygon.vertices.size(); j++) {
            Vector cur_vertex = clipped_polygon.vertices[j];
            Vector prev_vertex = clipped_polygon.vertices[(j > 0) ? (j - 1) : (clipped_polygon.vertices.size() - 1)];

            Vector intersection = intersect(prev_vertex, cur_vertex, point_u, point_v);

            if (inside(cur_vertex, point_u, point_v)) {
                if (!inside(prev_vertex, point_u, point_v)) {
                    out_polygon.vertices.push_back(intersection);
                }
                out_polygon.vertices.push_back(cur_vertex);
            } else if (inside(prev_vertex, point_u, point_v)) {
                out_polygon.vertices.push_back(intersection);
            }
        }
        clipped_polygon = out_polygon;
    }
    return clipped_polygon;
}

// Voronoi Parallel Linear Enumeration
Vector voronoi_intersect(const Vector &point_A, const Vector &point_B,
                         const Vector &point_Pi, const Vector &point_Pj)
{
    Vector mid_Pi_Pj = (point_Pi + point_Pj) * 0.5;
    double distance = dot(mid_Pi_Pj - point_A, point_Pi - point_Pj) / dot(point_B - point_A, point_Pi - point_Pj);
    if (!(distance >= 0 && distance <= 1)) return point_A;
    return point_A + (point_B - point_A) * distance;
}

bool voronoi_inside(const Vector &point_P, const Vector &point_Pi, const Vector &point_Pj) {
    Vector mid_Pi_Pj = (point_Pi + point_Pj) * 0.5;
    return dot(point_P - mid_Pi_Pj, point_Pj - point_Pi) < 0;
}

std::vector<Polygon> voronoi_ple(const std::vector<Vector> &sites) {   
    std::vector<Polygon> voronoi_cells;
    voronoi_cells.reserve(sites.size());

    for (int i = 0; i < sites.size(); i++) {
        const Vector &point_Pi = sites[i];

        // Largely shaped quadrilateral that contains all sites (unit square)
        Polygon clipped_polygon({
            Vector(0.0, 0.0),
            Vector(1.0, 0.0),
            Vector(1.0, 1.0),
            Vector(0.0, 1.0)
        });

        for (int j = 0; j < sites.size(); j++) {
            if (j == i) continue;
            const Vector &point_Pj = sites[j];

            Polygon out_polygon;

            for (int k = 0; k < clipped_polygon.vertices.size(); k++) {
                const Vector &cur_vertex  = clipped_polygon.vertices[k];
                const Vector &prev_vertex = clipped_polygon.vertices[(k > 0) ? (k - 1) : (clipped_polygon.vertices.size() - 1)];

                Vector intersection = voronoi_intersect(prev_vertex, cur_vertex, point_Pi, point_Pj);

                if (voronoi_inside(cur_vertex,  point_Pi, point_Pj)) {
                    if (!voronoi_inside(prev_vertex, point_Pi, point_Pj)) {
                        out_polygon.vertices.push_back(intersection);
                    }
                    out_polygon.vertices.push_back(cur_vertex);
                } else if (voronoi_inside(prev_vertex, point_Pi, point_Pj)) {
                    out_polygon.vertices.push_back(intersection);
                }
            }
            clipped_polygon = out_polygon;
        }
        voronoi_cells.push_back(clipped_polygon);
    }
    return voronoi_cells;
}

// Weighted Voronoi Parallel Linear Enumeration
Vector weighted_voronoi_intersect(const Vector &point_A, const Vector &point_B,
                         const Vector &point_Pi, const Vector &point_Pj,
                         const double &weight_i, const double &weight_j)
{
    Vector weighted_mid_Pi_Pj = (point_Pi + point_Pj) * 0.5 + ((weight_i - weight_j) / (2 * (dot(point_Pi - point_Pj, point_Pi - point_Pj)))) * (point_Pj - point_Pi);
    double distance = dot(weighted_mid_Pi_Pj - point_A, point_Pi - point_Pj) / dot(point_B - point_A, point_Pi - point_Pj);
    if (!(distance >= 0 && distance <= 1)) return point_A;
    return point_A + (point_B - point_A) * distance;
}

bool weighted_voronoi_inside(const Vector &point_P, 
                    const Vector &point_Pi, const Vector &point_Pj,
                    const double &weight_i, const double &weight_j) {
    Vector weighted_mid_Pi_Pj = (point_Pi + point_Pj) * 0.5 + ((weight_i - weight_j) / (2 * (dot(point_Pi - point_Pj, point_Pi - point_Pj)))) * (point_Pj - point_Pi);
    return dot(point_P - weighted_mid_Pi_Pj, point_Pj - point_Pi) < 0;
}

std::vector<Polygon> weighted_voronoi_ple(const std::vector<Vector> &sites, const std::vector<double> &weights) {   
    std::vector<Polygon> voronoi_cells;
    voronoi_cells.reserve(sites.size());

    for (int i = 0; i < sites.size(); i++) {
        const Vector &point_Pi = sites[i];

        // Largely shaped quadrilateral that contains all sites
        Polygon clipped_polygon({
            Vector(0.0, 0.0),
            Vector(1.0, 0.0),
            Vector(1.0, 1.0),
            Vector(0.0, 1.0)
        });

        for (int j = 0; j < sites.size(); j++) {
            if (j == i) continue;
            const Vector &point_Pj = sites[j];

            Polygon out_polygon;

            for (int k = 0; k < clipped_polygon.vertices.size(); k++) {
                const Vector &cur_vertex  = clipped_polygon.vertices[k];
                const Vector &prev_vertex = clipped_polygon.vertices[(k > 0) ? (k - 1) : (clipped_polygon.vertices.size() - 1)];

                Vector intersection = weighted_voronoi_intersect(prev_vertex, cur_vertex, point_Pi, point_Pj, weights[i], weights[j]);

                if (weighted_voronoi_inside(cur_vertex,  point_Pi, point_Pj, weights[i], weights[j])) {
                    if (!weighted_voronoi_inside(prev_vertex, point_Pi, point_Pj, weights[i], weights[j])) {
                        out_polygon.vertices.push_back(intersection);
                    }
                    out_polygon.vertices.push_back(cur_vertex);
                } else if (weighted_voronoi_inside(prev_vertex, point_Pi, point_Pj, weights[i], weights[j])) {
                    out_polygon.vertices.push_back(intersection);
                }
            }
            clipped_polygon = out_polygon;
        }
        voronoi_cells.push_back(clipped_polygon);
    }
    return voronoi_cells;
}

Polygon clip_by_disk(const Polygon& polygon, const Vector& disk_center, double disk_radius) {
    if (disk_radius <= 0 || polygon.vertices.empty()) 
        return polygon;
    
    // We consider a regular polygon with many vertices in order to simulate a circle
    int sides = 256;
    Polygon disk;
    for (int i = 0; i < sides; ++i) {
        double angle = 2 * M_PI * i / sides;
        disk.vertices.emplace_back(
            disk_center[0] + disk_radius * std::cos(angle),
            disk_center[1] + disk_radius * std::sin(angle)
        );
    }
    
    return polygon_clipping(polygon, disk);
}
