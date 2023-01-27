/*
 * Copyright 2020 Improbable Worlds Limited
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "benchmark_util.h"
#include "logging.h"
#include "phtree/phtree.h"
#include "phtree/phtree_multimap.h"
#include "phtree/phtree_multimap2.h"
#include "src/boost/boost_multimap.h"
#include "src/flann/ph-kdtree.h"
#include "src/tinspin/kdtree.h"
#include "src/tinspin/quadtree_point.h"
#include <benchmark/benchmark.h>
#include <random>

using namespace improbable;
using namespace improbable::phtree;
using namespace improbable::phtree::phbenchmark;

/*
 * Benchmark for k-nearest-neighbour queries in multi-map implementations.
 */
namespace {

const double GLOBAL_MAX = 10000;

enum Scenario {
    BOOST_RT,
    LSI,
    TREE_WITH_MAP,
    MULTI_MAP,
    MULTI_MAP_STD,
    PHTREE2,
    TS_KD,
    TS_QT,
    FLANN_KD
};

using payload_t = int64_t;
using payload2_t = uint32_t;

using TestPoint = PhPointD<3>;
using QueryBox = PhBoxD<3>;
using BucketType = std::set<payload_t>;

template <dimension_t DIM>
using CONVERTER = ConverterIEEE<DIM>;

template <Scenario SCENARIO, dimension_t DIM>
using TestMap = typename std::conditional_t<
    SCENARIO == BOOST_RT,
    b::PhTreeMultiMapD<DIM, payload_t>,
    typename std::conditional_t<
        //        SCENARIO == LSI,
        //        si::PhTreeMultiMapD<DIM, payload_t>,
        //        typename std::conditional_t<
        SCENARIO == TREE_WITH_MAP,
        PhTreeD<DIM, BucketType, CONVERTER<DIM>>,
        typename std::conditional_t<
            SCENARIO == MULTI_MAP,
            PhTreeMultiMap<DIM, payload_t, CONVERTER<DIM>, b_plus_tree_hash_set<payload_t>>,
            typename std::conditional_t<
                SCENARIO == PHTREE2,
                PhTreeMultiMap2D<DIM, payload_t>,
                typename std::conditional_t<
                    SCENARIO == MULTI_MAP_STD,
                    PhTreeMultiMap<DIM, payload_t, CONVERTER<DIM>, BucketType>,
                    typename std::conditional_t<
                        //                    SCENARIO == BB,
                        //                    bb::PhTreeMultiMapD<DIM, payload2_t>,
                        //                    typename std::conditional_t<
                        SCENARIO == TS_KD,
                        tinspin::KDTree<payload, double>,
                        typename std::conditional_t<
                            SCENARIO == TS_QT,
                            tinspin::QuadTree<payload>,
                            typename std::conditional_t<
                                SCENARIO == FLANN_KD,
                                flann::PhTreeMultiMap<DIM, size_t>,
                                void>>>>>>>>;

template <dimension_t DIM, Scenario SCENARIO>
class IndexBenchmark {
  public:
    IndexBenchmark(benchmark::State& state, int knn_result_size_);

    void Benchmark(benchmark::State& state);

  private:
    void SetupWorld(benchmark::State& state);
    void QueryWorld(benchmark::State& state, TestPoint& center);
    void CreateQuery(TestPoint& center);

    const TestGenerator data_type_;
    const size_t num_entities_;
    const size_t knn_result_size_;

    TestMap<SCENARIO, DIM> tree_;
    std::default_random_engine random_engine_;
    std::uniform_real_distribution<> cube_distribution_;
    std::vector<TestPoint> points_;
};

template <dimension_t DIM, Scenario SCENARIO>
IndexBenchmark<DIM, SCENARIO>::IndexBenchmark(benchmark::State& state, int knn_result_size)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, knn_result_size_(knn_result_size)
, random_engine_{1}
, cube_distribution_{0, GLOBAL_MAX}
, points_(state.range(0)) {
    logging::SetupDefaultLogging();
    SetupWorld(state);
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::Benchmark(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        TestPoint center;
        CreateQuery(center);
        state.ResumeTiming();

        QueryWorld(state, center);
    }
}

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<(SCENARIO != Scenario::TREE_WITH_MAP), int> = 0>
void InsertEntry(TestMap<SCENARIO, DIM>& tree, const TestPoint& point, const payload_t& data) {
    tree.emplace(point, data);
}

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<(SCENARIO == Scenario::TREE_WITH_MAP), int> = 0>
void InsertEntry(
    TestMap<Scenario::TREE_WITH_MAP, DIM>& tree, const TestPoint& point, const payload_t& data) {
    BucketType& bucket = tree.emplace(point).first;
    bucket.emplace(data);
}

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<(SCENARIO != Scenario::FLANN_KD), int> = 0>
size_t QueryAll(TestMap<SCENARIO, DIM>& tree, const TestPoint& center, const size_t k) {
    size_t n = 0;
    for (auto q = tree.begin_knn_query(k, center, DistanceEuclidean<3>()); q != tree.end(); ++q) {
        ++n;
    }
    return n;
}

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<(SCENARIO == Scenario::FLANN_KD), int> = 0>
size_t QueryAll(TestMap<SCENARIO, DIM>& tree, const TestPoint& center, const size_t k) {
    size_t n = 0;
    for (auto q = tree.begin_knn_query(k, center, DistanceEuclidean<3>()); q != tree.knn_end();
         ++q) {
        ++n;
    }
    return n;
}

struct CounterTreeWithMap {
    void operator()(const TestPoint&, const BucketType& value) {
        for (auto& x : value) {
            (void)x;
            n_ += 1;
        }
    }
    size_t n_;
};

struct CounterMultiMap {
    void operator()(const TestPoint&, const payload_t&) {
        n_ += 1;
    }
    size_t n_;
};

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    CreatePointData<DIM>(points_, data_type_, num_entities_, 0, GLOBAL_MAX);
    for (size_t i = 0; i < num_entities_; ++i) {
        InsertEntry<DIM, SCENARIO>(tree_, points_[i], (payload_t)i);
    }

    state.counters["query_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["result_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["avg_result_count"] = benchmark::Counter(0, benchmark::Counter::kAvgIterations);
    logging::info("World setup complete.");
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::QueryWorld(benchmark::State& state, TestPoint& center) {
    size_t n = QueryAll<DIM, SCENARIO>(tree_, center, knn_result_size_);

    state.counters["query_rate"] += 1;
    state.counters["result_rate"] += n;
    state.counters["avg_result_count"] += n;
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::CreateQuery(TestPoint& center) {
    for (dimension_t d = 0; d < DIM; ++d) {
        center[d] = cube_distribution_(random_engine_) * GLOBAL_MAX;
    }
}

}  // namespace

template <typename... Arguments>
void BoostRT(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::BOOST_RT> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// template <typename... Arguments>
// void LibSI(benchmark::State& state, Arguments&&... arguments) {
//     IndexBenchmark<3, Scenario::LSI> benchmark{state, arguments...};
//     benchmark.Benchmark(state);
// }
//
// template <typename... Arguments>
// void BBTree3D(benchmark::State& state, Arguments&&... arguments) {
//    IndexBenchmark<3, Scenario::BB> benchmark{state, arguments...};
//    benchmark.Benchmark(state);
//}

template <typename... Arguments>
void KDTree3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::TS_KD> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void Quadtree3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::TS_QT> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void FlannKD3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::FLANN_KD> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTree3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::TREE_WITH_MAP> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMultiMap3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::MULTI_MAP> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMultiMap2_3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::PHTREE2> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMultiMapStd3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::MULTI_MAP_STD> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_type, num_entities, query_result_size
// PhTree 3D CUBE

BENCHMARK_CAPTURE(PhTree3D, KNN_1, 1)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTree3D, KNN_10, 10)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

// Multimap 2.0
BENCHMARK_CAPTURE(PhTreeMultiMap2_3D, KNN_1, 1)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTreeMultiMap2_3D, KNN_10, 10)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

// FLANN KD-tree
BENCHMARK_CAPTURE(FlannKD3D, KNN_1, 1)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(FlannKD3D, KNN_10, 10)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

// KD-tree
BENCHMARK_CAPTURE(KDTree3D, KNN_1, 1)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(KDTree3D, KNN_10, 10)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

// Quadtree
BENCHMARK_CAPTURE(Quadtree3D, KNN_1, 1)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(Quadtree3D, KNN_10, 10)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 100 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

// BENCHMARK_CAPTURE(PhTree3D, KNN_CU_1_of_10K, TestGenerator::CUBE, 10000, 1)
//     ->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(PhTree3D, KNN_CU_1_of_1M, TestGenerator::CUBE, 1000000, 1)
//     ->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(PhTree3D, KNN_CU_10_of_10K, TestGenerator::CUBE, 10000, 10)
//     ->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(PhTree3D, KNN_CU_10_of_1M, TestGenerator::CUBE, 1000000, 10)
//     ->Unit(benchmark::kMillisecond);
//
//// index type, scenario name, data_type, num_entities, query_result_size
//// PhTree 3D CLUSTER
// BENCHMARK_CAPTURE(PhTree3D, KNN_CL_1_of_10K, TestGenerator::CLUSTER, 10000, 1)
//     ->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(PhTree3D, KNN_CL_1_of_1M, TestGenerator::CLUSTER, 1000000, 1)
//     ->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(PhTree3D, KNN_CL_10_of_10K, TestGenerator::CLUSTER, 10000, 10)
//     ->Unit(benchmark::kMillisecond);
//
// BENCHMARK_CAPTURE(PhTree3D, KNN_CL_10_of_1M, TestGenerator::CLUSTER, 1000000, 10)
//     ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
