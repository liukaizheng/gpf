# Rust to C++ Conversion Guide for GPF

This document describes the conversion of the GPF (Geometry Processing Functions) library from Rust to C++.

## Project Structure

### Rust Structure
```
src/
├── lib.rs              # Main library entry
├── math/
│   └── mod.rs
├── disjoint_set/
│   └── mod.rs
├── graphcut/
│   ├── mod.rs
│   └── push_relabel_fifo.rs
├── predicates/
│   ├── mod.rs
│   ├── expansion_number.rs
│   ├── generic_point.rs
│   ├── interval_number.rs
│   ├── less_than.rs
│   ├── orient2d.rs
│   ├── orient3d.rs
│   └── predicates.rs
├── triangle/
│   ├── mod.rs
│   ├── triangle.rs
│   └── tetrahedron.rs
├── mesh/
│   ├── mod.rs
│   ├── mesh.rs
│   ├── mesh_core_data.rs
│   ├── manifold_mesh.rs
│   ├── surface_mesh.rs
│   └── element/
│       ├── mod.rs
│       ├── vertex.rs
│       ├── edge.rs
│       ├── face.rs
│       ├── halfedge.rs
│       └── boundary_loop.rs
└── polygonlization/
    ├── mod.rs
    ├── bsp_complex.rs
    └── conforming_mesh.rs
```

### C++ Structure
```
include/gpf/          # Public headers
├── gpf.hpp           # Main header
├── math.hpp
├── disjoint_set.hpp
├── graphcut.hpp
├── predicates.hpp
├── triangle.hpp
├── mesh.hpp
└── polygonlization.hpp

src/                  # Implementation files
├── graphcut.cpp
├── predicates.cpp
├── triangle.cpp
├── mesh.cpp
└── polygonlization.cpp

tests/               # Unit tests
├── CMakeLists.txt
├── test_math.cpp
├── test_disjoint_set.cpp
├── test_graphcut.cpp
└── test_predicates.cpp
```

## Key Conversion Patterns

### 1. Module System

**Rust:**
```rust
pub mod math;
pub use math::*;
```

**C++:**
```cpp
// In header: include/gpf/math.hpp
namespace gpf {
namespace math {
    // declarations
}
}

// In main header: include/gpf/gpf.hpp
#include "gpf/math.hpp"
```

### 2. Ownership and Borrowing

**Rust:**
```rust
fn process(data: &[f64]) -> Vec<f64>
fn process_mut(data: &mut Vec<f64>)
fn process_owned(data: Vec<f64>)
```

**C++:**
```cpp
std::vector<double> process(const double* data, size_t n);
void process_mut(std::vector<double>& data);
void process_owned(std::vector<double> data);
```

### 3. Traits to Concepts/Templates

**Rust:**
```rust
trait GenericNum = Sized + Add<Output = Self> + Sub<Output = Self>;

fn compute<T: GenericNum>(a: T, b: T) -> T {
    a + b
}
```

**C++:**
```cpp
template<typename T>
concept GenericNum = requires(T a, T b) {
    { a + b } -> std::same_as<T>;
    { a - b } -> std::same_as<T>;
};

// Or with templates
template<typename T>
T compute(T a, T b) {
    return a + b;
}
```

### 4. Option and Result

**Rust:**
```rust
fn find(data: &[i32], target: i32) -> Option<usize>
fn parse(s: &str) -> Result<i32, ParseError>
```

**C++:**
```cpp
#include <optional>
#include <variant>

std::optional<size_t> find(const std::vector<int>& data, int target);

// For Result, use std::variant or std::expected (C++23)
std::variant<int, ParseError> parse(const std::string& s);
// Or exceptions
int parse(const std::string& s);  // throws ParseError
```

### 5. Iterators

**Rust:**
```rust
for item in vec.iter() {
    // process item
}

let result: Vec<_> = vec.iter().map(|x| x * 2).collect();
```

**C++:**
```cpp
for (const auto& item : vec) {
    // process item
}

std::vector<int> result;
std::transform(vec.begin(), vec.end(), std::back_inserter(result),
               [](int x) { return x * 2; });

// Or C++20 ranges
auto result = vec | std::views::transform([](int x) { return x * 2; })
                  | std::ranges::to<std::vector>();
```

### 6. Custom Allocators

**Rust (with Bumpalo):**
```rust
use bumpalo::Bump;
let arena = Bump::new();
let vec: Vec<f64, &Bump> = Vec::new_in(&arena);
```

**C++:**
```cpp
// Standard allocators
std::vector<double> vec;

// Custom allocator (if needed)
template<typename T>
class ArenaAllocator {
    // ... implementation
};

std::vector<double, ArenaAllocator<double>> vec;
```

### 7. HashMap

**Rust:**
```rust
use hashbrown::HashMap;
let mut map: HashMap<usize, Vec<usize>> = HashMap::new();
map.entry(key).or_insert(vec![]).push(value);
```

**C++:**
```cpp
#include <unordered_map>
std::unordered_map<size_t, std::vector<size_t>> map;
map[key].push_back(value);
```

### 8. Enums

**Rust:**
```rust
pub enum Orientation {
    Positive,
    Negative,
    Zero,
    Undefined,
}
```

**C++:**
```cpp
enum class Orientation {
    Positive,
    Negative,
    Zero,
    Undefined
};
```

### 9. Pattern Matching

**Rust:**
```rust
match orientation {
    Orientation::Positive => handle_positive(),
    Orientation::Negative => handle_negative(),
    _ => handle_other(),
}
```

**C++:**
```cpp
switch (orientation) {
    case Orientation::Positive:
        handle_positive();
        break;
    case Orientation::Negative:
        handle_negative();
        break;
    default:
        handle_other();
        break;
}
```

### 10. Inline Functions

**Rust:**
```rust
#[inline(always)]
fn fast_func(x: f64) -> f64 {
    x * x
}
```

**C++:**
```cpp
inline double fast_func(double x) {
    return x * x;
}

// Or force inline
[[gnu::always_inline]] inline double fast_func(double x) {
    return x * x;
}
```

## Build System Conversion

### Rust (Cargo.toml)
```toml
[package]
name = "gpf"
version = "0.1.0"
edition = "2021"

[dependencies]
bumpalo = "3.13.0"
hashbrown = "0.14.3"
itertools = "0.10.5"
rand = "0.8.5"
```

### C++ (CMakeLists.txt)
```cmake
cmake_minimum_required(VERSION 3.15)
project(gpf VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add library
add_library(gpf ${SOURCES})

# Add tests
enable_testing()
add_subdirectory(tests)
```

## Module-Specific Conversions

### Math Module
- **Status:** ✅ Complete
- **Changes:** Direct translation with template functions
- **Notes:** No major differences, straightforward conversion

### Disjoint Set
- **Status:** ✅ Complete
- **Changes:** 
  - `Vec<usize>` → `std::vector<size_t>`
  - `HashMap` → `std::unordered_map`
- **Notes:** Path compression implemented identically

### Graph Cut
- **Status:** ✅ Complete
- **Changes:**
  - Struct with public/private methods
  - Arrays for fixed-size queues
  - Vector growth instead of arena allocation
- **Notes:** Boykov-Kolmogorov algorithm preserved

### Predicates
- **Status:** ⚠️ Partial (basic predicates implemented)
- **Changes:**
  - Basic orient2d/orient3d/incircle/insphere implemented
  - Expansion arithmetic not yet ported
  - Filter/exact arithmetic cascade simplified
- **Todo:**
  - Port expansion number arithmetic
  - Implement interval arithmetic
  - Add robust predicates from Shewchuk

### Triangle
- **Status:** 🚧 Stub only
- **Todo:**
  - Port Delaunay triangulation
  - Port constrained Delaunay
  - Port 3D tetrahedralization
  - Implement incremental insertion

### Mesh
- **Status:** 🚧 Stub only
- **Todo:**
  - Complete half-edge data structure
  - Implement mesh iterators
  - Port mesh operations
  - Add manifold checking

### Polygonlization
- **Status:** 🚧 Stub only
- **Todo:**
  - Port BSP complex splitting
  - Implement conforming mesh generation
  - Add duplicate removal
  - Port mesh repair algorithm

## Performance Considerations

### Memory Management
- **Rust:** Uses Bumpalo arena allocators for hot paths
- **C++:** Uses standard allocators (can add custom allocators if needed)

### Optimization Flags
- **Rust:** `cargo build --release` uses `-O3` equivalent
- **C++:** Use `-O3 -march=native -DNDEBUG` for release builds

### Inlining
- **Rust:** `#[inline(always)]` for critical paths
- **C++:** `inline` and `[[gnu::always_inline]]` for critical paths

### SIMD
- **Rust:** `#![feature(portable_simd)]`
- **C++:** Consider using intrinsics or libraries like xsimd

## Testing

### Running Tests

**Rust:**
```bash
cargo test
cargo test --release
```

**C++:**
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
ctest -V
```

### Current Test Coverage
- ✅ Math operations
- ✅ Disjoint set union-find
- ✅ Graph cut max-flow
- ✅ Basic predicates (orient2d, orient3d, incircle, insphere)

## Next Steps

1. **Complete Predicates Module**
   - Implement expansion arithmetic
   - Add filter cascade for robust predicates
   - Port all geometric tests

2. **Triangle Module**
   - Port 2D Delaunay triangulation
   - Add constraint handling
   - Implement 3D tetrahedralization

3. **Mesh Module**
   - Complete half-edge implementation
   - Add mesh iterators (vertices, edges, faces)
   - Implement Euler operations

4. **Polygonlization Module**
   - Port BSP tree construction
   - Implement intersection splitting
   - Add mesh repair pipeline

5. **Optimization**
   - Profile and optimize hot paths
   - Consider custom allocators
   - Add SIMD where beneficial

## Known Limitations

1. **Expansion Arithmetic:** Not yet implemented - affects precision of predicates
2. **Arena Allocators:** Using standard allocators - may be slower for some operations
3. **Iterator Ergonomics:** C++ iterators less convenient than Rust iterators
4. **Memory Safety:** Manual memory management requires care
5. **Error Handling:** Using exceptions instead of Result types

## Resources

- Original Rust code in `src/` directory (preserved)
- C++ headers in `include/gpf/`
- C++ implementations in `src/`
- Tests in `tests/`
- Build with CMake

## License

Same as original Rust project.
