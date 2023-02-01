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
#include "src/bb-tree/bb_multimap.h"
#include "src/flann/ph-kdtree-single.h"
#include "phtree/phtree_multimap.h"
#include "phtree/phtree_multimap2.h"
#include "src/boost/boost_multimap.h"
#include "src/lsi/lsi_multimap.h"
#include "src/mcxme/mcxme_map.h"
#include "src/tinspin/kdtree.h"
#include "src/tinspin/quadtree_point.h"
#include <benchmark/benchmark.h>

using namespace improbable;
using namespace improbable::phtree;
using namespace improbable::phtree::phbenchmark;

namespace {

const double GLOBAL_MAX = 10000;

enum Scenario {
    BOOST_RT,
    LSI,
    MCXME,
    PHTREE,
    PHTREE2,
    TS_KD,
    TS_QT,
    BB,
    FLANN_KD_S,
};

using payload_t = int64_t;
using payload2_t = uint32_t;
using TestPointD = PhPointD<3>;
using TestPointF = std::vector<float>;
using TestPoint = TestPointD;

using BucketType = std::set<payload_t>;

template <dimension_t DIM>
using CONVERTER = ConverterIEEE<DIM>;

template <Scenario SCENARIO, dimension_t DIM>
using TestMap = typename std::conditional_t<
    SCENARIO == BOOST_RT,
    b::PhTreeMultiMapD<DIM, payload_t>,
    typename std::conditional_t<
        SCENARIO == LSI,
        si::PhTreeMultiMapD<DIM, payload_t>,
        typename std::conditional_t<
            SCENARIO == MCXME,
            mcxme::PhTreeD<DIM, payload_t>,
            typename std::conditional_t<
                SCENARIO == PHTREE,
                improbable::phtree::PhTreeMultiMapD<DIM, payload_t>,
                typename std::conditional_t<
                    SCENARIO == PHTREE2,
                    improbable::phtree::PhTreeMultiMap2D<DIM, payload_t>,
                    typename std::conditional_t<
                        SCENARIO == TS_KD,
                        tinspin::KDTree<TestPoint, payload_t>,
                        typename std::conditional_t<
                            SCENARIO == TS_QT,
                            tinspin::QuadTree<TestPoint, payload_t>,
                          typename std::conditional_t<
                              SCENARIO == BB,
                              bb::PhTreeMultiMapF<DIM, payload2_t>,
                              typename std::conditional_t<
                                  SCENARIO == FLANN_KD_S,
                                  flann::KDTreeSingle<DIM>,
                            void>>>>>>>>>;

/*
 * Benchmark for adding entries to the index.
 */
template <dimension_t DIM, Scenario S>
class IndexBenchmark {
    using Index = TestMap<S, DIM>;

  public:
    explicit IndexBenchmark(benchmark::State& state);
    void Benchmark(benchmark::State& state);

  private:
    void SetupWorld(benchmark::State& state);
    void Insert(benchmark::State& state, Index& tree);

    const TestGenerator data_type_;
    const size_t num_entities_;
  public:
    std::vector<TestPointD> points_d_;
    std::vector<TestPointF> points_f_;
    std::vector<std::uint32_t> values_;
};

template <dimension_t DIM, Scenario S>
IndexBenchmark<DIM, S>::IndexBenchmark(benchmark::State& state)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, points_d_(state.range(0)), points_f_(), values_() {
    logging::SetupDefaultLogging();
    SetupWorld(state);
}

template <dimension_t DIM, Scenario S>
void IndexBenchmark<DIM, S>::Benchmark(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        auto* tree = new Index();
        state.ResumeTiming();

        Insert(state, *tree);

        // we do this top avoid measuring deallocation
        state.PauseTiming();
        delete tree;
        state.ResumeTiming();
    }
}

template <dimension_t DIM, Scenario SCN>
void IndexBenchmark<DIM, SCN>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    CreatePointData<DIM>(points_d_, data_type_, num_entities_, 0, GLOBAL_MAX, 0.1);
    points_f_.reserve(num_entities_);
    values_.reserve(num_entities_);
    for (size_t i = 0; i < points_d_.size(); ++i) {
        values_.emplace_back(i);
        auto & v = points_f_.emplace_back(DIM);
        for (dimension_t d = 0; d < DIM; ++d) {
            v[d] = points_d_[i][d];
        }
    }

    state.counters["total_put_count"] = benchmark::Counter(0);
    state.counters["put_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    logging::info("World setup complete.");
}

//template <
//    dimension_t DIM,
//    Scenario SCN,
//    std::enable_if_t<(SCN != Scenario::BB && SCN != Scenario::FLANN_KD_S), int> = 0>
//void InsertEntries(TestMap<SCN, DIM>& tree, const std::vector<TestPoint>& points) {
//    for (size_t i = 0; i < points.size(); ++i) {
//        auto& p = points[i];
//        tree.emplace(p, (payload_t)i);
//    }
//}
//
//template <
//    dimension_t DIM,
//    Scenario SCN,
//    std::enable_if_t<(SCN == Scenario::BB || SCN == Scenario::FLANN_KD_S), int> = 0>
//void InsertEntries(TestMap<SCN, DIM>& tree, const std::vector<TestPoint>& points) {
//    tree.load(points);
//}
//
template <
    dimension_t DIM,
    Scenario SCN,
    std::enable_if_t<(SCN != Scenario::BB && SCN != Scenario::FLANN_KD_S), int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    const auto & points = data.points_d_;
    for (size_t i = 0; i < points.size(); ++i) {
        auto& p = points[i];
        tree.emplace(p, (payload_t)i);
    }
}

template <
    dimension_t DIM,
    Scenario SCN,
    std::enable_if_t<(SCN == Scenario::FLANN_KD_S), int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    tree.load(data.points_d_);
}

template <
    dimension_t DIM,
    Scenario SCN,
    std::enable_if_t<(SCN == Scenario::BB), int> = 0>
void InsertEntries(TestMap<SCN, DIM>& tree, const IndexBenchmark<DIM, SCN>& data) {
    tree.load(data.points_f_, data.values_);
}

template <dimension_t DIM, Scenario SCN>
void IndexBenchmark<DIM, SCN>::Insert(benchmark::State& state, Index& tree) {
    InsertEntries<DIM, SCN>(tree, *this);
    //    for (size_t i = 0; i < num_entities_; ++i) {
    //        auto& p = points_[i];
    //        tree.emplace(p, (payload_t)i);
    //    }

    state.counters["total_put_count"] += num_entities_;
    state.counters["put_rate"] += num_entities_;
}

}  // namespace

template <typename... Arguments>
void BoostRT(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, BOOST_RT> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void Lsi(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, LSI> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void Mcxme(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, MCXME> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMM(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, PHTREE> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void PhTreeMM2(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, PHTREE2> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void TinspinKDTree(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, TS_KD> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void TinspinQuadtree(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, TS_QT> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void BBTree(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, BB> benchmark{state};
    benchmark.Benchmark(state);
}

template <typename... Arguments>
void FlannKDS(benchmark::State& state, Arguments&&...) {
    IndexBenchmark<3, FLANN_KD_S> benchmark{state};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_generator, num_entities
BENCHMARK_CAPTURE(PhTreeMM, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 10}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(BBTree, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(FlannKDS, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(PhTreeMM2, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(TinspinKDTree, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(TinspinQuadtree, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(BoostRT, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(Mcxme, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(Lsi, INSERT, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
