# GPF C++ Quick Start Guide

## 5-Minute Setup

### Prerequisites
```bash
# Ubuntu/Debian
sudo apt install build-essential cmake

# macOS
brew install cmake

# Windows (use MSVC or MinGW)
# Download CMake from cmake.org
```

### Build in 3 Commands
```bash
mkdir build && cd build
cmake ..
make -j4
```

### Run Tests
```bash
ctest
```

Expected output:
```
Test project /path/to/build
    Start 1: test_math
1/4 Test #1: test_math ........................   Passed
    Start 2: test_disjoint_set
2/4 Test #2: test_disjoint_set ................   Passed
    Start 3: test_graphcut
3/4 Test #3: test_graphcut ....................   Passed
    Start 4: test_predicates
4/4 Test #4: test_predicates ..................   Passed

100% tests passed, 0 tests failed out of 4
```

### Run Example
```bash
./example
```

Expected output:
```
=== GPF C++ Library Demo ===

1. Math Operations
   Subtraction: [-3, -3, -3]
   Dot product: 32
   Cross product: [-3, 6, -3]

2. Disjoint Set (Union-Find)
   Initial groups: 6
   After merging pairs: 3 groups
   After connecting components: 2 groups

3. Graph Cut (Max Flow)
   Maximum flow: 37
   ...

4. Geometric Predicates
   Orient2D(p1, p2, p3): 0.5 (counter-clockwise)
   ...
```

## Use in Your Project

### Option 1: CMake Subproject

Add to your `CMakeLists.txt`:
```cmake
add_subdirectory(gpf)
target_link_libraries(your_target gpf)
```

### Option 2: Install System-Wide

```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
sudo make install
```

Then in your project:
```cmake
find_package(gpf REQUIRED)
target_link_libraries(your_target gpf)
```

### Option 3: Header-Only (for small modules)

Copy `include/gpf/` to your project and include:
```cpp
#include "gpf/math.hpp"
#include "gpf/disjoint_set.hpp"
```

## Basic Usage Examples

### Math Operations
```cpp
#include <gpf/math.hpp>

double a[] = {1, 2, 3};
double b[] = {4, 5, 6};

// Dot product
double d = gpf::math::dot(a, b, 3);  // 32

// Cross product
double c[3];
gpf::math::cross(a, b, c);  // [-3, 6, -3]

// Norm
double n = gpf::math::norm(a, 3);  // sqrt(14)
```

### Disjoint Set (Union-Find)
```cpp
#include <gpf/disjoint_set.hpp>

gpf::DisjointSet ds(10);  // 10 elements

// Union operations
ds.merge(0, 1);
ds.merge(2, 3);
ds.merge(1, 2);  // Now 0,1,2,3 are connected

// Get connected components
auto groups = ds.output();
std::cout << "Number of groups: " << ds.n_groups << "\n";
```

### Graph Cut (Max Flow)
```cpp
#include <gpf/graphcut.hpp>

std::vector<double> source_cap = {10, 5, 15};
std::vector<double> sink_cap = {8, 10, 7};

gpf::GraphCut gc(source_cap, sink_cap);

// Add edges (node_i, node_j, capacity, reverse_capacity)
gc.add_edge(0, 1, 10, 10);
gc.add_edge(1, 2, 5, 5);

// Compute max flow
double flow = gc.max_flow();
std::cout << "Max flow: " << flow << "\n";

// Check which nodes are in source/sink sets
for (size_t i = 0; i < gc.is_sink.size(); ++i) {
    if (gc.is_sink[i]) {
        std::cout << "Node " << i << " is in sink set\n";
    }
}
```

### Geometric Predicates
```cpp
#include <gpf/predicates.hpp>

// 2D orientation test
double p1[] = {0, 0};
double p2[] = {1, 0};
double p3[] = {0.5, 0.5};

double orient = gpf::predicates::orient2d(p1, p2, p3);
if (orient > 0) {
    std::cout << "Counter-clockwise\n";
} else if (orient < 0) {
    std::cout << "Clockwise\n";
} else {
    std::cout << "Collinear\n";
}

// 3D orientation test
double q1[] = {0, 0, 0};
double q2[] = {1, 0, 0};
double q3[] = {0, 1, 0};
double q4[] = {0, 0, 1};

double orient3 = gpf::predicates::orient3d(q1, q2, q3, q4);

// Incircle test
double result = gpf::predicates::incircle(p1, p2, p3, q4);
```

## What's Implemented

✅ **Fully Working:**
- Math operations (dot, cross, norm, sub)
- Disjoint Set / Union-Find
- Graph Cut (Boykov-Kolmogorov max-flow)
- Basic predicates (orient2d, orient3d, incircle, insphere)

🚧 **Stubs (need implementation):**
- Triangle (Delaunay triangulation)
- Mesh (half-edge mesh structures)
- Polygonlization (mesh repair)

## Common Issues

### Issue: CMake not found
```bash
# Ubuntu/Debian
sudo apt install cmake

# macOS
brew install cmake
```

### Issue: C++20 not supported
Ensure you have a modern compiler:
- GCC 10+
- Clang 11+
- MSVC 2019+

### Issue: Tests fail
```bash
# Clean rebuild
cd build
rm -rf *
cmake ..
make -j4
ctest -V  # Verbose output
```

## Performance Tips

### Release Build
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

### Optimizations
```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -march=native" ..
```

### Disable Tests
```bash
cmake -DBUILD_TESTS=OFF ..
```

### Disable Examples
```bash
cmake -DBUILD_EXAMPLES=OFF ..
```

## Documentation

- `README_CPP.md` - Full API documentation
- `CONVERSION_GUIDE.md` - Rust to C++ conversion patterns
- `CPP_PORT_SUMMARY.md` - Project status and statistics
- Headers in `include/gpf/` - API reference (well commented)
- Tests in `tests/` - Usage examples

## Next Steps

1. ✅ Build and run tests (you are here!)
2. 📖 Read `README_CPP.md` for full API docs
3. 💡 Check `example.cpp` for more examples
4. 🔧 Integrate into your project
5. 🚀 Start using GPF!

## Help & Support

- Check the documentation files listed above
- Look at test files for usage examples
- Compare with original Rust code in `src/` directories
- File issues if you find bugs

## License

Same as the original Rust project.

---

**Happy coding with GPF C++!** 🎉
