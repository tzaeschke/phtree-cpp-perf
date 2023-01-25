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
#include "src/lsi/lsi_multimap.h"
#include "src/robin_hood/robin_hood.h"
#include "src/tinspin/kdtree.h"
#include "src/tinspin/quadtree_point.h"
#include <benchmark/benchmark.h>
#include <random>

using namespace improbable;
using namespace improbable::phtree;
using namespace improbable::phtree::phbenchmark;

/*
 * Benchmark for querying entries in multi-map implementations.
 */
namespace {

const double GLOBAL_MAX = 10000;
const double BOX_LEN = 10;
const double QUERY_SHELL_RADIUS = BOX_LEN * 0.1;

const int DUMMY = 0;

enum Scenario {
    BOOST_RT,
    LSI,
    TREE_WITH_MAP,
    MULTI_MAP,
    MULTI_MAP_STD,
    MULTIPLY_MAP,
    PHTREE2,
    TS_KD,
    TS_QT
};

using payload_t = int64_t;

using TestPoint = PhPointD<3>;
using QueryBox = PhBoxD<3>;
using BucketType = std::set<payload_t>;
// using BucketType = std::unordered_set<payload_t>;
// using BucketType = robin_hood::unordered_set<payload_t>;

template <dimension_t DIM>
using CONVERTER = ConverterIEEE<DIM>;

template <dimension_t DIM>
using CONVERTER_MUL = ConverterMultiply<DIM, 1, 30>;

template <Scenario SCENARIO, dimension_t DIM>
using TestMap = typename std::conditional_t<
    SCENARIO == BOOST_RT,
    b::PhTreeMultiMapD<DIM, payload_t>,
    typename std::conditional_t<
        SCENARIO == LSI,
        si::PhTreeMultiMapD<DIM, payload_t>,
        typename std::conditional_t<
            SCENARIO == TREE_WITH_MAP,
            PhTreeD<DIM, BucketType, CONVERTER<DIM>>,
            typename std::conditional_t<
                SCENARIO == MULTI_MAP,
                PhTreeMultiMap<DIM, payload_t, CONVERTER<DIM>, b_plus_tree_hash_set<payload_t>>,
                typename std::conditional_t<
                    SCENARIO == MULTIPLY_MAP,
                    PhTreeMultiMapD<DIM, payload_t, CONVERTER_MUL<DIM>>,
                    typename std::conditional_t<
                        SCENARIO == PHTREE2,
                        PhTreeMultiMap2D<DIM, payload_t>,
                        typename std::conditional_t<
                            SCENARIO == MULTI_MAP_STD,
                            PhTreeMultiMap<DIM, payload_t, CONVERTER<DIM>, BucketType>,
                            typename std::conditional_t<
                                SCENARIO == TS_KD,
                                tinspin::KDTree<payload_t>,
                                typename std::conditional_t<
                                    SCENARIO == TS_QT,
                                    tinspin::QuadTree<payload_t>,
                                    void>>>>>>>>>;

template <dimension_t DIM, Scenario SCENARIO>
class IndexBenchmark {
  public:
    IndexBenchmark(benchmark::State& state, int);

    void Benchmark(benchmark::State& state);

  private:
    void SetupWorld(benchmark::State& state);
    void QueryWorld(benchmark::State& state, const QueryBox& query);
    void CreateQuery(QueryBox& query);

    const TestGenerator data_type_;
    const size_t num_entities_;

  public:
    TestMap<SCENARIO, DIM> tree_;
    std::default_random_engine random_engine_;
    std::uniform_int_distribution<> id_distribution_;
    std::vector<TestPoint> points_;
};

template <dimension_t DIM, Scenario SCENARIO>
IndexBenchmark<DIM, SCENARIO>::IndexBenchmark(benchmark::State& state, int)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, tree_{}
, random_engine_{1}
, id_distribution_{0, static_cast<int>(num_entities_)}
, points_(state.range(0)) {
    logging::SetupDefaultLogging();
    SetupWorld(state);
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::Benchmark(benchmark::State& state) {
    QueryBox query{};
    for (auto _ : state) {
        state.PauseTiming();
        CreateQuery(query);
        state.ResumeTiming();

        QueryWorld(state, query);
    }
}

template <dimension_t DIM, Scenario SCENARIO>
typename std::enable_if<SCENARIO == Scenario::TREE_WITH_MAP, void>::type InsertEntry(
    TestMap<Scenario::TREE_WITH_MAP, DIM>& tree, const TestPoint& point, const payload_t& data) {
    BucketType& bucket = tree.emplace(point).first;
    bucket.emplace(data);
}

template <dimension_t DIM, Scenario SCENARIO>
typename std::enable_if<SCENARIO != Scenario::TREE_WITH_MAP, void>::type InsertEntry(
    TestMap<SCENARIO, DIM>& tree, const TestPoint& point, const payload_t& data) {
    tree.emplace(point, data);
}

struct CounterTreeWithMap {
    void operator()(const TestPoint&, const BucketType& value) {
        for (auto it = value.begin(); it != value.end(); ++it) {
            n_ += 1;
        }
    }
    size_t n_;
};

struct CounterMultiMap {
    void operator()(const TestPoint&, const payload_t& value) {
        (void)value;
        n_ += 1;
    }
    size_t n_;
};

template <dimension_t DIM, Scenario SCENARIO>
typename std::enable_if<SCENARIO == Scenario::TREE_WITH_MAP, size_t>::type CountEntries(
    TestMap<Scenario::TREE_WITH_MAP, DIM>& tree, const QueryBox& query) {
    CounterTreeWithMap counter{0};
    tree.for_each(query, counter);
    return counter.n_;
}

template <dimension_t DIM, Scenario SCENARIO>
typename std::enable_if<SCENARIO != Scenario::TREE_WITH_MAP, size_t>::type CountEntries(
    TestMap<SCENARIO, DIM>& tree, const QueryBox& query) {
    CounterMultiMap counter{0};
    tree.for_each(query, counter);
    return counter.n_;
}

// template <dimension_t DIM, Scenario SCENARIO>
// size_t CountEntries(TestMap<Scenario::BOOST_RT, DIM>& tree, const QueryBox& query) {
//     CounterMultiMap counter{0};
//     tree.for_each(query.box, counter);
//     return counter.n_;
// }
//
// template <dimension_t DIM, Scenario SCENARIO>
// size_t CountEntries(TestMap<Scenario::MULTI_MAP, DIM>& tree, const QueryBox& query) {
//     CounterMultiMap counter{0};
//     tree.for_each(query.box, counter);
//     return counter.n_;
// }
//
// template <dimension_t DIM, Scenario SCENARIO>
// size_t CountEntries(TestMap<Scenario::MULTI_MAP_STD, DIM>& tree, const QueryBox& query) {
//     CounterMultiMap counter{0};
//     tree.for_each(query.box, counter);
//     return counter.n_;
// }

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    // create data with a point distance of BOX_LEN -> TODO we are reusing the argument -> hack!
    CreatePointData<DIM>(points_, data_type_, num_entities_, 0, GLOBAL_MAX, BOX_LEN);
    for (size_t i = 0; i < num_entities_; ++i) {
        InsertEntry<DIM, SCENARIO>(tree_, points_[i], (payload_t)i);
    }

    state.counters["query_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["result_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["avg_result_count"] = benchmark::Counter(0, benchmark::Counter::kAvgIterations);
    logging::info("World setup complete.");
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::QueryWorld(benchmark::State& state, const QueryBox& query) {
    size_t n = CountEntries<DIM, SCENARIO>(tree_, query);

    state.counters["query_rate"] += 1;
    state.counters["result_rate"] += n;
    state.counters["avg_result_count"] += n;
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::CreateQuery(QueryBox& query) {
    size_t id = id_distribution_(random_engine_);
    auto& point = points_[id];
    // Use BOX_LEN i.o. BOX_LEN/2 because we are querying for (center-)points, not boxes.
    double radius = BOX_LEN + QUERY_SHELL_RADIUS;
    for (dimension_t d = 0; d < DIM; ++d) {
        query.min()[d] = point[d] - radius;
        query.max()[d] = point[d] + radius;
    }
}

}  // namespace

template <typename... Arguments>
void BoostRT(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::BOOST_RT> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void LibSI(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::LSI> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void TinspinKDTree(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::TS_KD> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void TinspinQuadtree(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::TS_QT> benchmark{state, arguments...};
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

template <typename... Arguments>
void PhTreeMultiMapMultiply3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::MULTI_MAP> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_type, num_entities, avg_query_result_size
// PhTree
BENCHMARK_CAPTURE(PhTree3D, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(PhTreeMultiMap3D, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(PhTreeMultiMap2_3D, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(TinspinKDTree, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(TinspinQuadtree, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(PhTreeMultiMapStd3D, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap with multiply converter
BENCHMARK_CAPTURE(PhTreeMultiMapMultiply3D, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(BoostRT, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(LibSI, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
