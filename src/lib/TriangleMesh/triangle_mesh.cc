#include "triangle_mesh.h"

// BoundingBox class methods
void BoundingBox::compute(
    std::vector<Vector>::iterator& start, 
    std::vector<Vector>::iterator& end
) {
    if (start == end) return;
    Bmin = Vector(1.0 * INT_MAX, 1.0 * INT_MAX, 1.0 * INT_MAX);
    Bmax = Vector(1.0 * INT_MIN, 1.0 * INT_MIN, 1.0 * INT_MIN);
    for (std::vector<Vector>::iterator it = start; it != end; ++it) {
        Bmin.set_x(std::min(Bmin.get_x(), (*it).get_x()));
        Bmin.set_y(std::min(Bmin.get_y(), (*it).get_y()));
        Bmin.set_z(std::min(Bmin.get_z(), (*it).get_z()));

        Bmax.set_x(std::max(Bmax.get_x(), (*it).get_x()));
        Bmax.set_y(std::max(Bmax.get_y(), (*it).get_y()));
        Bmax.set_z(std::max(Bmax.get_z(), (*it).get_z()));
    }
}

void BoundingBox::compute(
    std::vector<TriangleIndices>::iterator& start,
    std::vector<TriangleIndices>::iterator& end,
    const std::vector<Vector>& vertices
) { 
    if (start == end) return;
    Bmin = Vector(1.0 * INT_MAX, 1.0 * INT_MAX, 1.0 * INT_MAX);
    Bmax = Vector(1.0 * INT_MIN, 1.0 * INT_MIN, 1.0 * INT_MIN);
    for (std::vector<TriangleIndices>::iterator it = start; it != end; ++it) {
        Vector triangle[3] = {
            vertices[(*it).vtxi],
            vertices[(*it).vtxj],
            vertices[(*it).vtxk]
        };
        for (int k = 0; k < 3; ++k) {
            Bmin.set_x(std::min(Bmin.get_x(), triangle[k].get_x()));
            Bmin.set_y(std::min(Bmin.get_y(), triangle[k].get_y()));
            Bmin.set_z(std::min(Bmin.get_z(), triangle[k].get_z()));
    
            Bmax.set_x(std::max(Bmax.get_x(), triangle[k].get_x()));
            Bmax.set_y(std::max(Bmax.get_y(), triangle[k].get_y()));
            Bmax.set_z(std::max(Bmax.get_z(), triangle[k].get_z()));
        }
    }
}

bool BoundingBox::intersect(Ray& ray, double& t_enter) {
    Vector vec_origin = ray.get_origin();
    Vector vec_unit_direction = ray.get_unit_direction();

    double txmin = (Bmin.get_x() - vec_origin.get_x()) / vec_unit_direction.get_x();
    double txmax = (Bmax.get_x() - vec_origin.get_x()) / vec_unit_direction.get_x();
    if (txmin > txmax) std::swap(txmin, txmax);

    double tymin = (Bmin.get_y() - vec_origin.get_y()) / vec_unit_direction.get_y();
    double tymax = (Bmax.get_y() - vec_origin.get_y()) / vec_unit_direction.get_y();
    if (tymin > tymax) std::swap(tymin, tymax);

    if ((txmin > tymax) || (tymin > txmax)) return false;

    txmin = std::max(txmin, tymin);
    txmax = std::min(txmax, tymax);

    double tzmin = (Bmin.get_z() - vec_origin.get_z()) / vec_unit_direction.get_z();
    double tzmax = (Bmax.get_z() - vec_origin.get_z()) / vec_unit_direction.get_z();
    if (tzmin > tzmax) std::swap(tzmin, tzmax);

    if ((txmin > tzmax) || (tzmin > txmax)) return false;
    
    t_enter = std::max(txmin, tzmin);

    return txmax >= std::max(t_enter, EPS);
}

// TriangleMesh class methods
TriangleMesh::TriangleMesh(
    Vector aux_vec_albedo, 
    bool aux_mirror,
    bool aux_transparent,
    bool aux_light_source,
    double aux_refraction_index,
    std::string aux_optimization
) {
    this->set_color(aux_vec_albedo);
    this->set_mirror(aux_mirror);
    this->set_transparent(aux_transparent);
    this->set_light_source(aux_light_source);
    this->set_refraction_index(aux_refraction_index);

    bounding_box_ready = false;
    bvh_ready = false;
    bvh = nullptr;

    optimization = aux_optimization;
}

void TriangleMesh::readOBJ(const char* obj) {
    char matfile[255];
    char grp[255];

    FILE* f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);

            } else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char* consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
                indices.push_back(t);
            } else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                } else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    } else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            } else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }

        }

    }
    fclose(f);
}

void TriangleMesh::compute_intersection(Ray& ray, const TriangleIndices& index, Intersection& intersection) {
    const Vector& A = vertices[index.vtxi];
    const Vector& B = vertices[index.vtxj];
    const Vector& C = vertices[index.vtxk];

    Vector vec_origin = ray.get_origin();
    Vector vec_unit_direction  = ray.get_unit_direction();

    Vector e1 = B - A;
    Vector e2 = C - A;

    double determinant = dot(vec_unit_direction, cross(e1, e2));

    double beta = dot(e2, cross(A - vec_origin, vec_unit_direction)) / determinant;
    double gamma = (-1) * dot(e1, cross(A - vec_origin, vec_unit_direction)) / determinant;
    double alpha = 1 - beta - gamma;
    double distance = dot(A - vec_origin, cross(e1, e2)) / determinant; 

    Vector vec_normal_A = normals[index.ni];
    Vector vec_normal_B = normals[index.nj];
    Vector vec_normal_C = normals[index.nk];
    Vector vec_normal = alpha * vec_normal_A + beta * vec_normal_B + gamma * vec_normal_C;
    vec_normal.normalize();

    bool cond_alpha    = (alpha >= 0) && (alpha <= 1);
    bool cond_beta     = (beta >= 0) && (beta <= 1);
    bool cond_gamma    = (gamma >= 0) && (gamma <= 1);
    bool cond_distance = (distance > EPS) && (distance < intersection.distance);

    if (cond_alpha && cond_beta && cond_gamma && cond_distance) {
        intersection.flag     = true;
        intersection.distance = distance;
        intersection.vec_point = vec_origin + vec_unit_direction * distance;
        intersection.vec_normal = vec_normal;
        intersection.vec_albedo = this->get_color();
    }
}

BVHNode* TriangleMesh::build_node(
    std::vector<TriangleIndices>::iterator start,
    std::vector<TriangleIndices>::iterator end
) {
    BVHNode* node = new BVHNode();

    node->bbox.compute(start, end, vertices);

    node->start_triangle = start;
    node->end_triangle   = end;

    Vector diag = node->bbox.Bmax - node->bbox.Bmin;
    Vector middle_diag = node->bbox.Bmin + diag * 0.5;

    int longest_axis = 0;
    if (diag.get_y() > diag.get_x()) longest_axis = 1;
    if (diag.get_z() > std::max(diag.get_x(), diag.get_y())) longest_axis = 2;

    auto pivot = start;
    for (std::vector<TriangleIndices>::iterator it = start; it != end; ++it) {
        Vector A = vertices[(*it).vtxi];
        Vector B = vertices[(*it).vtxj];
        Vector C = vertices[(*it).vtxk];
        Vector barycenter = (A + B + C) * (1.0 / 3.0);
        if (barycenter[longest_axis] < middle_diag[longest_axis]) {
            std::swap(*it, *pivot);
            ++pivot;
        }
    }

    if (pivot == start || pivot == end || std::distance(start, end) < 5) {
        node->left_bbox  = nullptr;
        node->right_bbox  = nullptr;
        return node;
    }

    node->left_bbox  = build_node(start, pivot);
    node->right_bbox  = build_node(pivot, end);
    return node;
}

void TriangleMesh::traverse_BVH(Ray& ray, BVHNode* node, Intersection& intersection) {
    double t_enter;
    if (!node || !node->bbox.intersect(ray, t_enter)) return;

    std::list<BVHNode*> nodes_to_visit;
    nodes_to_visit.push_front(node);

    double best_inter_distance = 1.0 * INT_MAX;

    while (!nodes_to_visit.empty()) {
        BVHNode* current_node = nodes_to_visit.back();
        nodes_to_visit.pop_back();
        double inter_distance = -1;

        if (current_node ->left_bbox) {
            if (current_node->left_bbox->bbox.intersect(ray, inter_distance)) {
                if (inter_distance < best_inter_distance) {
                    nodes_to_visit.push_back(current_node->left_bbox);
                }
            }

            if (current_node->right_bbox->bbox.intersect(ray, inter_distance)) {
                if (inter_distance < best_inter_distance) {
                    nodes_to_visit.push_back(current_node->right_bbox);
                }
            }
        }
        else {
            for (auto it = current_node->start_triangle; it != current_node->end_triangle; ++it) {
                compute_intersection(ray, *it, intersection);
            }
        }
    }
}

void TriangleMesh::build_BVH() {
    bvh = new BVHTree();
    bvh->root = build_node(indices.begin(), indices.end());
    bvh_ready = true;
}

void TriangleMesh::build_bbox() {
    std::vector<Vector>::iterator start = vertices.begin();
    std::vector<Vector>::iterator end = vertices.end();
    bounding_box.compute(start, end);
    bounding_box_ready = true;
}

Intersection TriangleMesh::intersected_by(Ray& ray) {
    Intersection intersection = Intersection();

    intersection.flag = false;
    intersection.distance = 1.0 * INT_MAX;

    if (optimization == "bvh" && bvh_ready) {
        traverse_BVH(ray, bvh->root, intersection);
    }
    else if (optimization == "bbox" || (optimization == "bvh" && !bvh_ready)) {
        double t_box;
        if (!bounding_box.intersect(ray, t_box)) {
            return intersection;
        }

        for (size_t i = 0; i < indices.size(); ++i) {
            compute_intersection(ray, indices[i], intersection);
        }
    }
    else if (optimization == "none" || (optimization == "bbox" && !bounding_box_ready)) {
        for (size_t i = 0; i < indices.size(); ++i) {
            compute_intersection(ray, indices[i], intersection);
        }
    }

    return intersection;
}
