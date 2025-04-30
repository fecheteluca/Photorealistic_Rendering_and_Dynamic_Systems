#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include <string>
#include <cstring>
#include <climits>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <atomic>
#include <list>
#include "vector.h"
#include "ray.h"
#include "geometry.h"

#define EPS 1e-6

class TriangleIndices {
    public:
        TriangleIndices(
            int vtxi = -1,
            int vtxj = -1,
            int vtxk = -1,
            int ni = -1,
            int nj = -1,
            int nk = -1,
            int uvi = -1,
            int uvj = -1,
            int uvk = -1,
            int group = -1
        ) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {}
    
        int vtxi, vtxj, vtxk;
        int uvi, uvj, uvk;
        int ni, nj, nk;
        int group;
    };

class BoundingBox {
public:
    Vector Bmin;
    Vector Bmax;

    void compute(
        std::vector<Vector>::iterator& start,
        std::vector<Vector>::iterator& end
    );
    void compute(
        std::vector<TriangleIndices>::iterator& start,
        std::vector<TriangleIndices>::iterator& end,
        const std::vector<Vector>& vertices
    );
    bool intersect(Ray& ray, double& t_enter);
};

class BVHNode {
public:
    BoundingBox bbox;
    BVHNode* left_bbox = nullptr;
    BVHNode* right_bbox = nullptr;

    std::vector<TriangleIndices>::iterator start_triangle;
    std::vector<TriangleIndices>::iterator end_triangle;
};

class BVHTree {
public:
    BVHNode* root = nullptr;
};

class TriangleMesh : public Geometry {
public:
    explicit TriangleMesh(
        Vector aux_albedo = Vector(),
        bool aux_mirror = false,
        bool aux_transparent = false,
        bool aux_light_source = false,
        double aux_refraction_index = 1.0,
        std::string aux_optimization = "bvh"
    );
    ~TriangleMesh() = default;

    void readOBJ(const char* obj);

    void compute_intersection(Ray& ray, const TriangleIndices& index, Intersection& intersection);

    BVHNode* build_node(std::vector<TriangleIndices>::iterator start,
                       std::vector<TriangleIndices>::iterator end);
    void traverse_BVH(Ray& ray, BVHNode* node, Intersection& best);
    void build_BVH();
    void build_bbox();

    virtual Intersection intersected_by(Ray& ray) override;

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

private:
    BoundingBox bounding_box;
    bool bounding_box_ready;

    BVHTree* bvh;
    bool bvh_ready;
    
    // optimization = "none" - this is the naive traversal method
    // optimization = "bbox" - this is the optimization using Bounding Boxes
    // optimization = "bvh"  - this is the optimization using BVH Trees
    std::string optimization;
};

#endif // TRIANGLEMESH_H
