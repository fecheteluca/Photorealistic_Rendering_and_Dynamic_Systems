#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include <string>
#include <cstring>
#include <climits>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "vector.h"
#include "ray.h"
#include "geometry.h"

#define EPS 1e-6

class BoundingBox {
public:
    Vector Bmin;
    Vector Bmax;

    void compute(std::vector<Vector>& verts);
    bool intersect(Ray& ray, double& t_enter);
};

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

class TriangleMesh : public Geometry {
public:
    explicit TriangleMesh(
        Vector aux_albedo = Vector(),
        bool aux_mirror = false,
        bool aux_transparent = false,
        bool aux_light_source = false,
        double aux_refraction_index = 1.0
    );
    ~TriangleMesh() = default;

    void readOBJ(const char* obj);

    virtual Intersection intersected_by(Ray& ray) override;

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;

private:
    BoundingBox bounding_box;
    bool bounding_box_ready;
};

#endif // TRIANGLEMESH_H
