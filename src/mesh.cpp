#include "gpf/mesh.hpp"
#include <algorithm>

namespace gpf {
namespace mesh {

SurfaceMesh::SurfaceMesh() {}

size_t SurfaceMesh::add_vertex(double x, double y, double z) {
    Vertex v;
    v.position[0] = x;
    v.position[1] = y;
    v.position[2] = z;
    v.halfedge = vertices_.size();
    vertices_.push_back(v);
    return vertices_.size() - 1;
}

size_t SurfaceMesh::add_face(const std::vector<size_t>& verts) {
    Face f;
    f.halfedge = halfedges_.size();
    f.vertices = verts;
    faces_.push_back(f);
    return faces_.size() - 1;
}

size_t SurfaceMesh::add_edge(size_t v1, size_t v2) {
    (void)v1; (void)v2;
    Edge e;
    e.halfedge = halfedges_.size();
    edges_.push_back(e);
    return edges_.size() - 1;
}

void SurfaceMesh::remove_vertex(size_t vid) {
    (void)vid;
}

void SurfaceMesh::remove_face(size_t fid) {
    (void)fid;
}

void SurfaceMesh::remove_edge(size_t eid) {
    (void)eid;
}

bool SurfaceMesh::is_manifold() const {
    return true;
}

bool SurfaceMesh::is_closed() const {
    return true;
}

ManifoldMesh::ManifoldMesh() : SurfaceMesh() {}

bool ManifoldMesh::check_manifold() {
    return true;
}

void ManifoldMesh::make_manifold() {
}

} // namespace mesh
} // namespace gpf
