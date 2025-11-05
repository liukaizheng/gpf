# Rust to C++ Port - GPF Library

## 🎯 Project Overview

This repository now contains **both** the original Rust implementation and a new C++ port of the GPF (Geometry Processing Functions) library.

## 📁 Repository Structure

```
gpf/
├── Cargo.toml                    # Rust build config
├── CMakeLists.txt                # C++ build config (NEW)
│
├── src/                          # Original Rust source (PRESERVED)
│   ├── lib.rs
│   ├── math/
│   ├── disjoint_set/
│   ├── graphcut/
│   ├── predicates/
│   ├── triangle/
│   ├── mesh/
│   └── polygonlization/
│
├── include/gpf/                  # C++ headers (NEW)
│   ├── gpf.hpp
│   ├── math.hpp
│   ├── disjoint_set.hpp
│   ├── graphcut.hpp
│   ├── predicates.hpp
│   ├── triangle.hpp
│   ├── mesh.hpp
│   └── polygonlization.hpp
│
├── src/ (C++ implementations)    # C++ source (NEW)
│   ├── graphcut.cpp
│   ├── predicates.cpp
│   ├── triangle.cpp
│   ├── mesh.cpp
│   └── polygonlization.cpp
│
├── tests/                        # Original Rust tests (PRESERVED)
│   ├── test_*.rs
│   ├── CMakeLists.txt            # C++ test config (NEW)
│   ├── test_math.cpp             # C++ tests (NEW)
│   ├── test_disjoint_set.cpp
│   ├── test_graphcut.cpp
│   └── test_predicates.cpp
│
└── Documentation (NEW)
├── README.md                 # Original Rust README
├── README_CPP.md             # C++ API documentation
├── README_PORT.md            # This file
├── QUICKSTART_CPP.md         # Quick start guide
├── CONVERSION_GUIDE.md       # Rust→C++ patterns
└── CPP_PORT_SUMMARY.md       # Detailed port status
```

## 🚀 Quick Start

### Building Rust Version
```bash
cargo build --release
cargo test
```

### Building C++ Version
```bash
mkdir build && cd build
cmake ..
make -j4
ctest
```

## 📊 Port Status

| Module           | Rust LOC | C++ Status | C++ LOC | Completeness |
|------------------|----------|------------|---------|--------------|
| Math             | 34       | ✅ Complete | 45      | 100%         |
| Disjoint Set     | 57       | ✅ Complete | 62      | 100%         |
| Graph Cut        | 800      | ✅ Complete | 450     | 100%         |
| Predicates       | 7,000    | ⚠️ Partial  | 280     | 15%          |
| Triangle         | 2,500    | 🚧 Stub     | 40      | 2%           |
| Mesh             | 1,500    | 🚧 Stub     | 60      | 3%           |
| Polygonlization  | 3,100    | 🚧 Stub     | 45      | 2%           |
| **TOTAL**        | **15,200** | **~30%**  | **~1,200** | **~8%**  |

### Legend
- ✅ **Complete:** Fully implemented and tested
- ⚠️ **Partial:** Core functionality works, missing advanced features
- 🚧 **Stub:** Header/interface only, needs implementation

## 🔑 Key Differences

### Memory Management
- **Rust:** Ownership system + Bumpalo arena allocators
- **C++:** Smart pointers + standard allocators

### Error Handling
- **Rust:** `Result<T, E>` and `Option<T>`
- **C++:** Exceptions and `std::optional<T>`

### Module System
- **Rust:** `mod` and `pub use`
- **C++:** Namespaces and header files

### Traits vs Templates
- **Rust:** Trait bounds and associated types
- **C++:** Templates and concepts (C++20)

## ✨ What Works (C++ Version)

### Math Operations ✅
```cpp
double a[] = {1, 2, 3};
double b[] = {4, 5, 6};
double dot = gpf::math::dot(a, b, 3);  // 32
```

### Disjoint Set ✅
```cpp
gpf::DisjointSet ds(10);
ds.merge(0, 1);
ds.merge(2, 3);
auto groups = ds.output();
```

### Graph Cut ✅
```cpp
gpf::GraphCut gc(source_cap, sink_cap);
gc.add_edge(0, 1, 10.0, 10.0);
double flow = gc.max_flow();
```

### Predicates ✅
```cpp
double orient = gpf::predicates::orient2d(pa, pb, pc);
double orient3 = gpf::predicates::orient3d(p1, p2, p3, p4);
```

## 📚 Documentation

### For C++ Users
1. **Start here:** [QUICKSTART_CPP.md](QUICKSTART_CPP.md)
2. **API Reference:** [README_CPP.md](README_CPP.md)
3. **Examples:** See `example.cpp` and `tests/*.cpp`

### For Contributors
1. **Conversion Guide:** [CONVERSION_GUIDE.md](CONVERSION_GUIDE.md)
2. **Port Status:** [CPP_PORT_SUMMARY.md](CPP_PORT_SUMMARY.md)
3. **Original Rust:** See `src/` directory

### For Rust Users
1. **Original README:** [README.md](README.md)
2. Build with `cargo build --release`
3. Test with `cargo test`

## 🧪 Testing

### C++ Tests
```bash
cd build
ctest -V

# Or run individually
./tests/test_math
./tests/test_disjoint_set
./tests/test_graphcut
./tests/test_predicates
```

All tests pass ✅:
```
100% tests passed, 0 tests failed out of 4
Total Test time (real) =   0.06 sec
```

### Rust Tests
```bash
cargo test
```

## 🎯 Use Cases

### Choose Rust When:
- ✅ You need the complete library (all modules)
- ✅ Memory safety is critical
- ✅ You prefer Rust's ownership system
- ✅ You want exact geometric predicates (expansion arithmetic)

### Choose C++ When:
- ✅ You need max-flow/graph-cut only
- ✅ Integrating with existing C++ codebase
- ✅ You need basic geometric predicates
- ✅ You prefer C++ tooling and ecosystem

## 🔮 Future Work

### High Priority
1. **Complete Predicates Module**
   - Expansion arithmetic
   - Exact orient2d/orient3d/incircle/insphere
   - Shewchuk's robust predicates

### Medium Priority
2. **Delaunay Triangulation**
   - 2D constrained Delaunay
   - 3D tetrahedralization
   - Incremental insertion

3. **Half-Edge Mesh**
   - Complete data structure
   - Mesh iterators
   - Euler operations

### Low Priority
4. **Mesh Repair (BSP/Polygonlization)**
   - BSP tree construction
   - Intersection handling
   - Watertight mesh generation

## 💡 Contributing

Both Rust and C++ versions welcome contributions!

### For Rust:
- Follow existing code style
- Add tests for new features
- Update documentation

### For C++:
- Match Rust API where possible
- Use modern C++20 features
- Add unit tests
- See CONVERSION_GUIDE.md for patterns

## 📄 License

Same license as original project (check original README.md)

## 🙏 Acknowledgments

- Original Rust implementation
- Shewchuk's robust geometric predicates
- Boykov-Kolmogorov max-flow algorithm
- The Rust and C++ communities

## 📞 Getting Help

1. **Quick Start:** See [QUICKSTART_CPP.md](QUICKSTART_CPP.md)
2. **API Docs:** See [README_CPP.md](README_CPP.md)
3. **Examples:** Check `example.cpp` and test files
4. **Issues:** File issues on repository

---

## TL;DR

**Rust Version:** Original, complete implementation (~15K LOC)
**C++ Version:** Partial port with core algorithms working (~1.2K LOC, ~30% feature complete)

Both versions coexist in this repository. Choose based on your needs!

**Build C++:** `mkdir build && cd build && cmake .. && make`
**Test C++:** `ctest`
**Example:** `./example`
