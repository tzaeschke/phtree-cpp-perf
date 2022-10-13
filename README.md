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
* CLUSTER dataset: A cube with clusters of 100 entries each. Clusters follow a Gauss distribution with standard deviation 0.0001.



## Benchmarks

Benchmarks filenames are more or less descriptive.
* `_mm_` indicates "multi-map" tests. Most spatial indexes are **multimaps** but the PH-tree works actually better as **map** (with no key duplicates)
* `hd_` indicates a test with > 3 dimensions, usually 6D, 10D and 20D
* `_box_` indicates a test with boxes (otherwise it is points)
* `_d_` indicates a test with `double` coordinates (otherwise it is integer)

## Output

Output it encoded. The first field is the index name, e.g.
1) Index name (and configuration), e.g.
   * `BoostRT` for boost R-Tree
   * `LSI` for libSpatialIndex R-Tree
   * `MCXME` for MCXME PH-tree
   * Others: PH-tree variants
2) Benchmark type, e.g.
   * `WQ` for window query
3) Number of entries in the index, usually powers of 10.
4) Dataset, e.g.
   * `4` CUBE dataset (4 letters)
   * `7` CLUSTER dataset (7 letters)
5) Optional: additional fields, e.g.
   * For window queries: expected average number of results per query window
   
The following example shows the result of a Boost R-Tree with a window query on an index with 1000 entries using the CLUSTER dataset:
```
BoostRT/WQ/1000/7                     0.000 ms        0.000 ms      2684182          2.59242 3.77716M/s  9.79197M/s
```




# Usage
The benchmarks can be run with bazel or cmake.

E.g. bazel:
```
bazel run //benchmark:insert_mm_box_d_benchmark --config=benchmark -- --benchmark_counters_tabular=true
```


# Dependencies

This project contains libraries form other projects:
* robin_hood 3.11.5, licensed under MIT: https://github.com/martinus/robin-hood-hashing

This projects depends on the following other projects:
* Boost (geometry), licensed under BSL-1.0: https://github.com/boostorg/geometry
* libSpatialIndex, licensed under MIT: https://github.com/libspatialindex/libspatialindex
* MCXME PH-tree, licensed under Apache 2.0: https://github.com/mcxme/phtree
* spdlog, licensed under MIT: https://github.com/gabime/spdlog
* GoogleTest, licensed under BSD 3-clause: https://github.com/google/googletest
* Bazel, licensed under Apache 2.0: https://github.com/bazelbuild/bazel
* CMake, licensed under BSD 3-clause: https://cmake.org


