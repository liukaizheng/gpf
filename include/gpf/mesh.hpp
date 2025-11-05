#pragma once

#include <vector>
#include <optional>
#include <memory>
#include <unordered_map>

namespace gpf {
namespace mesh {

struct HalfEdge {
    size_t vertex;
    size_t face;
    size_t next;
    size_t prev;
    size_t twin;
};

struct Vertex {
    double position[3];
    size_t halfedge;
};

struct Face {
    size_t halfedge;
    std::vector<size_t> vertices;
};

struct Edge {
    size_t halfedge;
};

class SurfaceMesh {
public:
    SurfaceMesh();
    
    size_t add_vertex(double x, double y, double z);
    size_t add_face(const std::vector<size_t>& vertices);
    size_t add_edge(size_t v1, size_t v2);
    
    void remove_vertex(size_t vid);
    void remove_face(size_t fid);
    void remove_edge(size_t eid);
    
    const std::vector<Vertex>& vertices() const { return vertices_; }
    const std::vector<Face>& faces() const { return faces_; }
    const std::vector<Edge>& edges() const { return edges_; }
    const std::vector<HalfEdge>& halfedges() const { return halfedges_; }
    
    size_t n_vertices() const { return vertices_.size(); }
    size_t n_faces() const { return faces_.size(); }
    size_t n_edges() const { return edges_.size(); }
    
    bool is_manifold() const;
    bool is_closed() const;
    
private:
    std::vector<Vertex> vertices_;
    std::vector<Face> faces_;
    std::vector<Edge> edges_;
    std::vector<HalfEdge> halfedges_;
};

class ManifoldMesh : public SurfaceMesh {
public:
    ManifoldMesh();
    
    bool check_manifold();
    void make_manifold();
};

} // namespace mesh
} // namespace gpf
