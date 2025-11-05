#pragma once

#include <vector>
#include <tuple>

namespace gpf {
namespace polygonlization {

std::tuple<std::vector<double>, std::vector<size_t>> 
make_mesh_for_triangles(const std::vector<double>& points,
                        const std::vector<size_t>& triangles,
                        const std::vector<size_t>& tri_in_shells);

class BSPComplex {
public:
    BSPComplex(const std::vector<double>& points, const std::vector<size_t>& triangles);
    
    void split_intersecting_faces();
    void remove_duplicates();
    std::tuple<std::vector<double>, std::vector<size_t>> extract_mesh();
    
private:
    std::vector<double> points_;
    std::vector<size_t> triangles_;
};

class ConformingMesh {
public:
    ConformingMesh(const std::vector<double>& points, const std::vector<size_t>& triangles);
    
    void make_conforming();
    std::tuple<std::vector<double>, std::vector<size_t>> get_mesh();
    
private:
    std::vector<double> points_;
    std::vector<size_t> triangles_;
};

} // namespace polygonlization
} // namespace gpf
