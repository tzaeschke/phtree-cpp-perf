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
#include <benchmark/benchmark.h>
#include <random>

using namespace improbable;
using namespace improbable::phtree;
using namespace improbable::phtree::phbenchmark;

/*
 * Benchmark for querying entries in multi-map implementations.
 */
namespace {
// TODO - check optimised query
// TODO - check libSI
const double GLOBAL_MAX = 10000;
const double BOX_LEN = 10;
const double QUERY_SHELL_RADIUS = BOX_LEN * 0.1;

const int DUMMY = 0;

enum Scenario {
    BOOST_RT,
    LSI_RT,
    TREE_SET,
    PHTREE_MM,
    MULTIPLY_MAP,
    PHTREE2,
    TS_KD,
    TS_QT,
    FLANN_KD_S
};

using payload_t = int64_t;

using TestPoint = PhPointD<3>;
using QueryBox = PhBoxD<3>;
using BucketType = std::unordered_set<payload_t>;

template <Scenario SCENARIO, dimension_t DIM>
using CONVERTER = ConverterBoxIEEE<DIM>;

template <Scenario SCENARIO, dimension_t DIM>
using CONVERTER_MUL = ConverterBoxMultiply<DIM, 1, 30>;

template <Scenario SCENARIO, dimension_t DIM>
using TestMap = typename std::conditional_t<
    SCENARIO == BOOST_RT,
    b::PhTreeMultiMapBoxD<DIM, payload_t>,
    typename std::conditional_t<
        SCENARIO == LSI_RT,
        si::PhTreeMultiMapBoxD<DIM, payload_t>,
        typename std::conditional_t<
            SCENARIO == TREE_SET,
            PhTreeBoxD<DIM, BucketType, CONVERTER<SCENARIO, DIM>>,
            typename std::conditional_t<
                SCENARIO == MULTIPLY_MAP,
                PhTreeMultiMapBoxD<DIM, BucketType, CONVERTER_MUL<SCENARIO, DIM>>,
                typename std::conditional_t<
                    SCENARIO == PHTREE2,
                    PhTreeMultiMap2BoxD<DIM, payload_t>,
                    PhTreeMultiMapBoxD<DIM, payload_t, CONVERTER<SCENARIO, DIM>>>>>>>;

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
    std::vector<PhBoxD<DIM>> boxes_;
};

template <dimension_t DIM, Scenario SCENARIO>
IndexBenchmark<DIM, SCENARIO>::IndexBenchmark(benchmark::State& state, int)
: data_type_{static_cast<TestGenerator>(state.range(1))}
, num_entities_(state.range(0))
, tree_{}
, random_engine_{1}
, id_distribution_{0, static_cast<int>(num_entities_)}
, boxes_(num_entities_) {
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

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<(SCENARIO != Scenario::TREE_SET), int> = 0>
void InsertEntry(TestMap<SCENARIO, DIM>& tree, const PhBoxD<DIM>& point, const payload_t& data) {
    tree.emplace(point, data);
}

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<(SCENARIO == Scenario::TREE_SET), int> = 0>
void InsertEntry(
    TestMap<Scenario::TREE_SET, DIM>& tree, const PhBoxD<DIM>& point, const payload_t& data) {
    BucketType& bucket = tree.emplace(point).first;
    bucket.emplace(data);
}

// bool CheckPosition(const PhBoxD<3>& box, const QueryBox& query) {
//     //   const auto& box = entity;
//     bool result = true;
//     for (int d = 0; d < 3; ++d) {
//         result = query.min()[d] <= box.max()[0] && query.max()[d] >= box.min()[d];
//     }
//     return result;
// }

struct CounterTreeWithMap {
    void operator()(const PhBoxD<3>&, const BucketType& value) {
        for (auto& x : value) {
            (void)x;
            n_ += 1;  // CheckPosition(boxes_[x], box_); // TODO
        }
    }
    // std::vector<PhBoxD<3>> boxes_;
    const QueryBox& box_;
    size_t n_;
};

struct CounterMultiMap {
    void operator()(const PhBoxD<3>&, const payload_t& value) {
        (void)value;
        n_ += 1;  // CheckPosition(boxes_[value], box_); // TODO
    }
    // std::vector<PhBoxD<3>> boxes_;
    const QueryBox& box_;
    size_t n_;
};

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<SCENARIO == Scenario::TREE_SET, int> = 0>
size_t CountEntries(TestMap<Scenario::TREE_SET, DIM>& tree, const QueryBox& query) {
    CounterTreeWithMap counter{query, 0};
    tree.for_each(query, counter);
    return counter.n_;
}

template <
    dimension_t DIM,
    Scenario SCENARIO,
    std::enable_if_t<SCENARIO != Scenario::TREE_SET, int> = 0>
size_t CountEntries(TestMap<SCENARIO, DIM>& tree, const QueryBox& query) {
    CounterMultiMap counter{query, 0};
    tree.for_each(query, counter);
    return counter.n_;
}

template <dimension_t DIM, Scenario SCENARIO>
void IndexBenchmark<DIM, SCENARIO>::SetupWorld(benchmark::State& state) {
    logging::info("Setting up world with {} entities and {} dimensions.", num_entities_, DIM);
    // create data with about 10% duplicate coordinates
    CreateBoxData<DIM>(boxes_, data_type_, num_entities_, 0, GLOBAL_MAX, BOX_LEN);
    for (size_t i = 0; i < num_entities_; ++i) {
        InsertEntry<DIM, SCENARIO>(tree_, boxes_[i], (payload_t)i);
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
    auto& box = boxes_[id];
    for (dimension_t d = 0; d < DIM; ++d) {
        query.min()[d] = box.min()[d] - QUERY_SHELL_RADIUS;
        query.max()[d] = box.max()[d] + QUERY_SHELL_RADIUS;
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
void PhTreeMultiMapMultiply3D(benchmark::State& state, Arguments&&... arguments) {
    IndexBenchmark<3, Scenario::PHTREE_MM> benchmark{state, arguments...};
    benchmark.Benchmark(state);
}

// index type, scenario name, data_type, num_entities, avg_query_result_size
// PhTree
BENCHMARK_CAPTURE(PhTreeSet, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(PhTreeMM, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::WEB, TestGenerator::WEB}})
    ->Unit(benchmark::kMillisecond);

// PhTreeMultiMap
BENCHMARK_CAPTURE(PhTreeMM2, WEB, DUMMY)
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

// LSI
BENCHMARK_CAPTURE(LsiRT, WEB, DUMMY)
    ->RangeMultiplier(10)
    ->Ranges({{1000, 1000 * 1000}, {TestGenerator::CUBE, TestGenerator::CLUSTER}})
    ->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
