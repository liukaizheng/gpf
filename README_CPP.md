# GPF C++ - Geometry Processing Functions

This is a C++ port of the GPF (Geometry Processing Functions) library, originally written in Rust. The library implements various geometric processing algorithms for robust solid modeling.

## Overview

GPF provides a collection of geometric processing functions with a focus on performance and correctness. The library includes:

- **Math Module**: Basic vector operations (dot, cross, norm, subtraction)
- **Disjoint Set**: Union-find data structure for connectivity analysis
- **Graph Cut**: Boykov-Kolmogorov max-flow/min-cut algorithm
- **Predicates**: Robust geometric predicates (orient2d, orient3d, incircle, insphere)
- **Triangle**: 2D Delaunay triangulation and 3D tetrahedralization
- **Mesh**: Half-edge mesh data structure for surface and manifold meshes
- **Polygonlization**: Mesh repair and watertight solid construction

## Building

### Requirements

- CMake 3.15 or higher
- C++20 compatible compiler (GCC 10+, Clang 11+, MSVC 2019+)

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

### Build with tests

```bash
mkdir build
cd build
cmake -DBUILD_TESTS=ON ..
make
ctest
```

## Usage

Include the main header:

```cpp
#include <gpf/gpf.hpp>
```

### Example: Math Operations

```cpp
#include <gpf/math.hpp>

double a[] = {1.0, 2.0, 3.0};
double b[] = {4.0, 5.0, 6.0};
double c[3];

// Subtraction
gpf::math::sub(a, b, c, 3);

// Norm
double norm = gpf::math::norm(a, 3);

// Dot product
double dot = gpf::math::dot(a, b, 3);

// Cross product
double cross[3];
gpf::math::cross(a, b, cross);
```

### Example: Disjoint Set

```cpp
#include <gpf/disjoint_set.hpp>

gpf::DisjointSet ds(5);  // 5 elements
ds.merge(0, 1);
ds.merge(2, 3);
auto groups = ds.output();  // Get connected components
```

### Example: Graph Cut

```cpp
#include <gpf/graphcut.hpp>

std::vector<double> source_cap = {10.0, 5.0, 15.0};
std::vector<double> sink_cap = {8.0, 10.0, 7.0};

gpf::GraphCut gc(source_cap, sink_cap);
gc.add_edge(0, 1, 10.0, 10.0);
gc.add_edge(1, 2, 5.0, 5.0);

double max_flow = gc.max_flow();
```

### Example: Geometric Predicates

```cpp
#include <gpf/predicates.hpp>

double pa[] = {0.0, 0.0};
double pb[] = {1.0, 0.0};
double pc[] = {0.0, 1.0};

// 2D orientation test
double orient = gpf::predicates::orient2d(pa, pb, pc);

// 3D orientation test
double p3d[] = {0.0, 0.0, 1.0};
double orient3 = gpf::predicates::orient3d(pa, pb, pc, p3d);

// Incircle test
double result = gpf::predicates::incircle(pa, pb, pc, pd);
```

## Module Status

This is a port from Rust to C++. The following modules have been implemented:

- ✅ Math: Complete
- ✅ Disjoint Set: Complete
- ✅ Graph Cut: Complete
- ✅ Predicates: Basic implementation (orient2d, orient3d, incircle, insphere)
- 🚧 Triangle: Stub implementation (requires full porting)
- 🚧 Mesh: Stub implementation (requires full porting)
- 🚧 Polygonlization: Stub implementation (requires full porting)

## Differences from Rust Version

1. **Memory Management**: C++ uses standard library containers instead of Bumpalo arena allocators
2. **Error Handling**: Uses exceptions and std::optional instead of Result types
3. **Traits**: Replaced with templates and abstract base classes where appropriate
4. **Ownership**: Uses smart pointers (unique_ptr, shared_ptr) for resource management

## Performance

The C++ version aims to match the performance characteristics of the Rust version:

- Inline functions for hot paths
- Efficient memory allocation patterns
- Cache-friendly data structures
- SIMD-friendly algorithms (where applicable)

## Testing

Run tests with:

```bash
cd build
ctest -V
```

Individual test executables:
- `test_math`: Math operations
- `test_disjoint_set`: Union-find
- `test_graphcut`: Max-flow algorithm
- `test_predicates`: Geometric predicates

## License

Same as the original Rust version.

## Contributing

This is a work-in-progress port. The following areas need completion:

1. Full implementation of Triangle module (Delaunay triangulation, tetrahedralization)
2. Full implementation of Mesh module (half-edge operations, iterators)
3. Full implementation of Polygonlization module (BSP complex, conforming mesh)
4. Expansion arithmetic for exact predicates
5. Comprehensive test suite

## Original Project

This is a C++ port of the Rust project available at the original repository.
The Rust version can be found in the same directory structure.
