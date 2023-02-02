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
#include "src/bb-tree/bb_multimap.h"
#include "src/boost/boost_multimap.h"
#include "src/flann/ph-kdtree-single.h"
#include "src/lsi/lsi_multimap.h"
#include "src/mcxme/mcxme_map.h"
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
    LSI_RT,
    TREE_SET,
    PHTREE_MM,
    PHTREE_MM_STD,
    MULTIPLY_MAP,
    PHTREE2,
    BB,
    TS_KD,
    TS_QT,
    FLANN_KD_S,
    MCXME,
};

using payload_t = int64_t;
using payload2_t = uint32_t;

using TestPointD = PhPointD<3>;
using TestPointF = std::vector<float>;
using TestPoint = TestPointD;
using QueryBox = PhBoxD<3>;
using BucketType = std::set<payload_t>;
// using BucketType = std::unordered_set<payload_t>;
// using BucketType = robin_hood::unordered_set<payload_t>;

template <dimension_t DIM>
using CONVERTER = ConverterIEEE<DIM>;

template <dimension_t DIM>
using CONVERTER_MUL = ConverterMultiply<DIM, 1, 30>;

template <Scenario SCENARIO, dimension_t DIM>
using TestMap = typename std::conditional_t < SCENARIO == BOOST_RT,
      b::PhTreeMultiMapD<DIM, payload_t>,
      typename std::conditional_t<
          SCENARIO == LSI_RT,
          si::PhTreeMultiMapD<DIM, payload_t>,
          typename std::conditional_t<
              SCENARIO == TREE_SET,
              PhTreeD<DIM, BucketType, CONVERTER<DIM>>,
              typename std::conditional_t<
                  SCENARIO == PHTREE_MM,
                  PhTreeMultiMap<DIM, payload_t, CONVERTER<DIM>, b_plus_tree_hash_set<payload_t>>,
                  typename std::conditional_t<
                      SCENARIO == MULTIPLY_MAP,
                      PhTreeMultiMapD<DIM, payload_t, CONVERTER_MUL<DIM>>,
                      typename std::conditional_t<
                          SCENARIO == PHTREE2,
                          PhTreeMultiMap2D<DIM, payload_t>,
                          typename std::conditional_t<
                              SCENARIO == PHTREE_MM_STD,
                              PhTreeMultiMap<DIM, payload_t, CONVERTER<DIM>, BucketType>,
                              typename std::conditional_t<
                                  SCENARIO == TS_KD,
                                  tinspin::KDTree<TestPoint, payload_t>,
                                  typename std::conditional_t<
                                      SCENARIO == TS_QT,
                                      tinspin::QuadTree<TestPoint, payload_t>,
                                      typename std::conditional_t<
                                          SCENARIO == FLANN_KD_S,
                                          flann::KDTreeSingle<DIM, size_t>,
                                          typename std::conditional_t<
                                              SCENARIO == BB,
                                              bb::PhTreeMultiMap<DIM, payload2_t>,
                                              void>>>>>>>>>>>;

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

  public:
    std::vector<TestPointD> points_d_;
    std::vector<TestPointF> points_f_;
    std::vector<std::uint32_t> values_;
};

template <dimension_t DIM, Scenario SCENARIO>
IndexBenchmark<DIM, SCENARIO>::IndexBenchmark(benchmark::State& state, int)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, tree_{}
, random_engine_{1}
, id_distribution_{0, static_cast<int>(num_entities_)}
, points_d_(state.range(0))
, points_f_()
, values_() {
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

//template <dimension_t DIM, Scenario SCENARIO>
//typename std::enable_if<SCENARIO == Scenario::TREE_SET, void>::type InsertEntry(
//    TestMap<Scenario::TREE_SET, DIM>& tree, const TestPoint& point, const payload_t& data) {
//    BucketType& bucket = tree.emplace(point).first;
//    bucket.emplace(data);
//}
//
//template <dimension_t DIM, Scenario SCENARIO>
//typename std::enable_if<SCENARIO != Scenario::TREE_SET, void>::type InsertEntry(
//    TestMap<SCENARIO, DIM>& tree, const TestPoint& point, const payload_t& data) {
//    tree.emplace(point, data);
//}

template <
    dimension_t DIM,
    Scenario SCN,
    std::enable_if_t<
        (SCN != Scenario::BB && SCN != Scenario::FLANN_KD_S && SCN != Scenario::TREE_SET),
        int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    const auto& points = data.points_d_;
    for (size_t i = 0; i < points.size(); ++i) {
        auto& p = points[i];
        tree.emplace(p, (payload_t)i);
    }
}

template <dimension_t DIM, Scenario SCN, std::enable_if_t<(SCN == Scenario::FLANN_KD_S), int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    tree.load(data.points_d_);
}

template <dimension_t DIM, Scenario SCN, std::enable_if_t<(SCN == Scenario::BB), int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    tree.load(data.points_f_, data.values_);
}

template <dimension_t DIM, Scenario SCN, std::enable_if_t<(SCN == Scenario::TREE_SET), int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    const auto& points = data.points_d_;
    for (size_t i = 0; i < points.size(); ++i) {
        auto& point = points[i];
        BucketType& bucket = tree.emplace(point).first;
        bucket.emplace((payload_t)i);
    }
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
typename std::enable_if<SCENARIO == Scenario::TREE_SET, size_t>::type CountEntries(
    TestMap<Scenario::TREE_SET, DIM>& tree, const TestPoint& min, const TestPoint& max) {
    CounterTreeWithMap counter{0};
    tree.for_each({min, max}, counter);
    return counter.n_;
}

template <dimension_t DIM, Scenario SCENARIO>
typename std::enable_if<SCENARIO != Scenario::TREE_SET, size_t>::type CountEntries(
    TestMap<SCENARIO, DIM>& tree, const TestPoint& min, const TestPoint& max) {
    CounterMultiMap counter{0};
    tree.for_each({min, max}, counter);
    return counter.n_;
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    // create data with a point distance of BOX_LEN -> TODO we are reusing the argument -> hack!
    CreatePointData<DIM>(points_d_, data_type_, num_entities_, 0, GLOBAL_MAX, BOX_LEN);
    points_f_.reserve(num_entities_);
    values_.reserve(num_entities_);
    for (size_t i = 0; i < points_d_.size(); ++i) {
        values_.emplace_back(i);
        auto& v = points_f_.emplace_back(DIM);
        for (dimension_t d = 0; d < DIM; ++d) {
            v[d] = points_d_[i][d];
        }
    }

    InsertEntries<DIM, SCENARIO>(tree_, *this);

    state.counters["query_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    //    state.counters["result_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    state.counters["avg_result_count"] = benchmark::Counter(0, benchmark::Counter::kAvgIterations);
    logging::info("World setup complete.");
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::QueryWorld(benchmark::State& state, const QueryBox& query) {
    size_t n = CountEntries<DIM, SCENARIO>(tree_, query.min(), query.max());

    state.counters["query_rate"] += 1;
    //    state.counters["result_rate"] += n;
    state.counters["avg_result_count"] += n;
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::CreateQuery(QueryBox& query) {
    size_t id = id_distribution_(random_engine_);
    auto& point = points_d_[id];
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
void LsiRT(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::LSI_RT> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void BBTree(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::BB> benchmark{state, arguments...};
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
void FlannKDS(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::FLANN_KD_S> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeSet(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::TREE_SET> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMM(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::PHTREE_MM> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMM2(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::PHTREE2> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMMStdSet(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::PHTREE_MM_STD> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMultiMapMultiply3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::PHTREE_MM> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_type, num_entities, avg_query_result_size
// PhTreeMultiMap
BENCHMARK_CAPTURE(PhTreeMM, WQ, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap2
BENCHMARK_CAPTURE(PhTreeMM2, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// KD-tree
BENCHMARK_CAPTURE(TinspinKDTree, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// Quadtree
BENCHMARK_CAPTURE(TinspinQuadtree, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(BoostRT, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// Flann KD-tree
BENCHMARK_CAPTURE(FlannKDS, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// BB-tree
BENCHMARK_CAPTURE(BBTree, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(LsiRT, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap with multiply converter
BENCHMARK_CAPTURE(PhTreeMultiMapMultiply3D, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
