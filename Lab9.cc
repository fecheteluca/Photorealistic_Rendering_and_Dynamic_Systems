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
#include <cmath>
#include <map>
#include <set>

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
        bool operator==(const Vector& a);
    
    private:
        double data[3];
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


class TriangleMesh {
public:
    explicit TriangleMesh();

    void readOBJ(const char* obj);

    std::vector<Vector> tutteEmbedding(int nbIter) const;

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
};

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

std::vector<Vector> TriangleMesh::tutteEmbedding(int nbIter) const {
    using Edge = std::pair<int,int>;
    auto makeEdge = [&](int a,int b) {
        return a < b ? Edge(a,b) : Edge(b,a);
    };

    std::map<Edge,int> edgeCount;
    for (auto &t : indices) {
        int a = t.vtxi, b = t.vtxj, c = t.vtxk;
        edgeCount[makeEdge(a,b)]++;
        edgeCount[makeEdge(b,c)]++;
        edgeCount[makeEdge(c,a)]++;
    }

    std::vector<Edge> boundaryEdges;
    for (auto &ec : edgeCount) {
        if (ec.second == 1)
            boundaryEdges.push_back(ec.first);
    }

    std::unordered_map<int, std::vector<int>> bAdj;
    for (auto &e : boundaryEdges) {
        bAdj[e.first].push_back(e.second);
        bAdj[e.second].push_back(e.first);
    }

    std::vector<int> boundary;
    boundary.reserve(boundaryEdges.size());
    int start = boundaryEdges.front().first;
    boundary.push_back(start);
    int prev = start;
    int cur  = bAdj[start][0];        
    while (cur != start) {
        boundary.push_back(cur);
        auto &nbrs = bAdj[cur];
        int next = (nbrs[0] == prev ? nbrs[1] : nbrs[0]);
        prev = cur;
        cur  = next;
    }

    int B = int(boundary.size());
    double s = 0.0;
    for (int i = 0; i < B; ++i) {
        int bi = boundary[i];
        int bj = boundary[(i+1)%B];
        s += (vertices[bi] - vertices[bj]).norm();
    }

    std::vector<Vector> pos = vertices;  
    for (auto &v : pos) v.set_z(0);

    double cs = 0.0;
    for (int i = 0; i < B; ++i) {
        int bi = boundary[i];
        int bj = boundary[(i+1)%B];
        double θ = 2.0 * M_PI * (cs / s);
        pos[bi][0] = std::cos(θ);
        pos[bi][1] = std::sin(θ);
        cs += (vertices[bi] - vertices[bj]).norm();
    }

    int N = vertices.size();
    std::vector<std::set<int>> tmpNb(N);
    for (auto &t : indices) {
        int a = t.vtxi, b = t.vtxj, c = t.vtxk;
        tmpNb[a].insert(b); tmpNb[a].insert(c);
        tmpNb[b].insert(a); tmpNb[b].insert(c);
        tmpNb[c].insert(a); tmpNb[c].insert(b);
    }
    std::vector<std::vector<int>> neighbors(N);
    for (int i = 0; i < N; ++i) {
        neighbors[i].assign(tmpNb[i].begin(), tmpNb[i].end());
    }

    std::vector<char> isBoundary(N, 0);
    for (int b : boundary) isBoundary[b] = 1;

    for (int iter = 0; iter < nbIter; ++iter) {
        std::vector<Vector> nextPos = pos;  
        for (int i = 0; i < N; ++i) {
            if (isBoundary[i]) continue;
            auto &nbr = neighbors[i];
            Vector sum(0,0,0);
            for (int j : nbr) sum = sum + pos[j];
            nextPos[i] = sum / double(nbr.size());
        }
        pos.swap(nextPos);
    }

    return pos;
}

int main() {
    TriangleMesh *goethe_mesh = new TriangleMesh();
    goethe_mesh->readOBJ("../../objects/goethe/goethe.obj");

    std::vector<Vector> goethe_embedded = goethe_mesh->tutteEmbedding(100);

    return 0;
}
