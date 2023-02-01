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
    EMPLACE,
    SQUARE_BR,
    TS_KD,
    TS_QT,
};

using payload_t = int64_t;
using TestPoint = PhPointD<3>;

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
                            void>>>>>>>;

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
    std::vector<TestPoint> points_;
};

template <dimension_t DIM, Scenario S>
IndexBenchmark<DIM, S>::IndexBenchmark(benchmark::State& state)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, points_(state.range(0)) {
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

template <dimension_t DIM, Scenario S>
void IndexBenchmark<DIM, S>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    CreatePointData<DIM>(points_, data_type_, num_entities_, 0, GLOBAL_MAX, 0.1);

    state.counters["total_put_count"] = benchmark::Counter(0);
    state.counters["put_rate"] = benchmark::Counter(0, benchmark::Counter::kIsRate);
    logging::info("World setup complete.");
}

template <dimension_t DIM, Scenario S>
void IndexBenchmark<DIM, S>::Insert(benchmark::State& state, Index& tree) {
    for (size_t i = 0; i < num_entities_; ++i) {
        auto& p = points_[i];
        tree.emplace(p, (payload_t)i);
    }

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

// index type, scenario name, data_generator, num_entities
BENCHMARK_CAPTURE(PhTreeMM, INSERT, 0)
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

BENCHMARK_CAPTURE(BoostRT, BOOST, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(Mcxme, MCXME, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_CAPTURE(Lsi, LSI, 0)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1 * 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
