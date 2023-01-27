# Spatial Index Benchmark

C++ Spatial index performance comparison.

This project compares performance of the [PH-tree](http://phtree.org) with various other spatial indexes.

Performance is compared with:

* Boost Geometry R-Tree: https://github.com/boostorg/geometry
* libSpatialIndex R-Tree: https://github.com/libspatialindex/libspatialindex
* MCXME PH-tree: https://github.com/mcxme/phtree

Operations that are compared:

* window query
* insert entry

Other operations show currently only PH-tree performance without comparison

* update position
* k-nearest-neighbor
* erase entry
* find entry
* ...

Comparisons use the following two datasets with varying size. Data ist mostly 3D. Datasets can be points or boxes.

* CUBE dataset: A cube with uniform randomly distributed points or boxes.
* CLUSTER dataset: A cube with clusters of 100 entries each. Clusters follow a Gauss distribution with standard
  deviation 0.0001.

## Benchmarks

Benchmarks filenames are more or less descriptive.

* `_mm_` indicates "multi-map" tests (otherwise it is **map**). Most spatial indexes are **multimaps** but the PH-tree
  works actually better as **map** (with no key duplicates)
* `hd_` indicates a test with > 3 dimensions, usually 6D, 10D and 20D
* `_box_` indicates a test with boxes (otherwise it is points)
* `_d_` indicates a test with `double` coordinates (otherwise it is `std::int64_t`)

## Output

Output it encoded. The fields have usually the following meaning:

1) Index name (and configuration), e.g.
    * `BoostRT` for boost R-Tree
    * `LSI` for libSpatialIndex R-Tree
    * `MCXME` for MCXME PH-tree
    * Others: PH-tree variants
2) Benchmark type, e.g.
    * `WQ` for window query
3) Number of entries in the index, usually powers of 10.
4) Dataset, e.g.
    * `4` for CUBE dataset (4 letters)
    * `7` for CLUSTER dataset (7 letters)
5) Optional: additional fields, e.g.
    * For window queries: expected average number of results per query window

The following example shows the result of a **Boost R-Tree** with a **window query** on an index with **1000** entries
using the **CLUSTER** dataset:

```
BoostRT/WQ/1000/7                     0.000 ms        0.000 ms      2684182          2.59242 3.77716M/s  9.79197M/s
```

## Caveats

Measuring performance for different APIs is not always easy.
This results in some necessary mapping from the PH-tree benchamarks to the index's APIs.

The Boost R-tree appeared to be the strongest contender of the three "foreign" indices.
It's mapping has been optimized considerable and should have only minimal impact.

Some notable issues:

* Boost RT is used with `std::array` coordinates instead of it's native coordinates. This appears to cause
  a small but measurable slowdown.
* liSpatialIndex uses a dirty hack for queries: results are buffered in a member field of the wrapper class.
* For `insert()`/`emplace()` the coordinate conversion for libSpatialIndex and MCXME may cost performance.
* There must be a bug: when using libSpatialIndex with WEB (box or points) it only finds 1 result on average where it
  should find 3 (as all other indexes do). This probelm does not occur with CUBE/CLUSTER or in tests.

# Usage

The benchmarks can be run with bazel or cmake.

E.g. bazel:

```
bazel run //benchmark:insert_mm_box_d_benchmark --config=benchmark -- --benchmark_counters_tabular=true
```

# Dependencies

This project contains libraries form other projects:

* MCXME PH-tree, licensed under Apache 2.0: https://github.com/mcxme/phtree
* PH-tree C++ 1.4.0, licensed under Apache 2.0: https://github.com/tzaeschke/phtree
  and https://github.com/improbable-eng/phtree-cpp
* robin_hood 3.11.5, licensed under MIT: https://github.com/martinus/robin-hood-hashing
* BB-tree, no license, Authors: Stefan Sprenger, Patrick Schäfer and Ulf Leser, Contact: sprengsz (at) informatik (dot)
  hu-berlin (dot) de: https://www2.informatik.hu-berlin.de/~sprengsz/bb-tree/
* PCL 1.13.0, licensed under BSD, Point Cloud Library (PCL): www.pointclouds.org
* FLANN 1.9.2, licensed under BSD, Fast Library for Approximate Nearest Neighbors: https://github.com/flann-lib/flann

This projects depends on the following other projects:

* Boost (geometry), licensed under BSL-1.0: https://github.com/boostorg/geometry
* libSpatialIndex, licensed under MIT: https://github.com/libspatialindex/libspatialindex
* spdlog, licensed under MIT: https://github.com/gabime/spdlog
* GoogleTest, licensed under BSD 3-clause: https://github.com/google/googletest
* Bazel, licensed under Apache 2.0: https://github.com/bazelbuild/bazel
* CMake, licensed under BSD 3-clause: https://cmake.org
* rules_boost, licensed under Apache 2.0: https://github.com/nelhage/rules_boost
* rules_pcl, licensed under Apache 2.0: https://github.com/kgreenek/rules_pcl

