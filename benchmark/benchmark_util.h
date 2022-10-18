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

#ifndef PHTREE_BENCHMARK_UTIL_H
#define PHTREE_BENCHMARK_UTIL_H

#include "phtree/common/common.h"
#include <random>
#include <vector>

namespace improbable::phtree::phbenchmark {

using namespace improbable::phtree;

namespace {
template <dimension_t DIM>
auto CreateDataCUBE = [](auto& points,
                         size_t num_entities,
                         std::uint32_t seed,
                         double world_mininum,
                         double world_maximum,
                         auto set_coordinate) {
    std::default_random_engine random_engine{seed};
    std::uniform_real_distribution<> distribution(world_mininum, world_maximum);
    for (size_t i = 0; i < num_entities; ++i) {
        auto& p = points[i];
        for (dimension_t d = 0; d < DIM; ++d) {
            set_coordinate(p, d, distribution(random_engine));
        }
    }
};

template <dimension_t DIM>
auto CreateDataCLUSTER = [](auto& points,
                            size_t num_entities,
                            std::uint32_t seed,
                            double world_mininum,
                            double world_maximum,
                            auto set_coordinate) {
    std::default_random_engine random_engine{seed};
    std::uniform_real_distribution<> distribution(world_mininum, world_maximum);
    // SIGMA = 0.0001
    std::normal_distribution<> gauss_distribution(0, 0.0001);
    const int NUM_PT_PER_CLUSTER = 100;
    // 1000 points per cluster, minimum is 1 cluster.
    size_t num_cluster = std::max(1, (int)(num_entities / NUM_PT_PER_CLUSTER));
    const double world_length = world_maximum - world_mininum;

    // loop over clusters
    PhPointD<DIM> cp;  // center point of cluster
    size_t id = 0;
    for (size_t c = 0; c < num_cluster; ++c) {
        for (dimension_t d = 0; d < DIM; ++d) {
            cp[d] = distribution(random_engine);
        }
        for (size_t i = 0; i < NUM_PT_PER_CLUSTER; ++i) {
            auto& p = points[id++];
            // int ii = (c * N_C + i) * DIM;
            for (dimension_t d = 0; d < DIM; ++d) {
                // double x = (R.nextGaussian() - 0.5) * GAUSS_SIGMA;  // confine to small rectangle
                double x = gauss_distribution(random_engine);
                x *= world_length + world_mininum;  // stretch if domain>1.0
                x += cp[d];                         // offset of cluster
                set_coordinate(p, d, x);
            }
        }
    }
};

template <dimension_t DIM>
using NodeD = std::array<double, DIM>;
template <dimension_t DIM>
double dist(const NodeD<DIM>& n1, const NodeD<DIM>& n2) {
    double sum = 0;
    for (dimension_t d = 0; d < DIM; ++d) {
        sum += (n1[d] - n2[d]) * (n1[d] - n2[d]);
    }
    return sqrt(sum);
}

template <dimension_t DIM>
void add(NodeD<DIM>& in_out, const NodeD<DIM>& diff) {
    for (dimension_t d = 0; d < DIM; ++d) {
        in_out[d] += diff[d];
    }
}

template <dimension_t DIM>
void sub(NodeD<DIM>& in_out, const NodeD<DIM>& diff) {
    for (dimension_t d = 0; d < DIM; ++d) {
        in_out[d] -= diff[d];
    }
}

template <dimension_t DIM>
void scale(NodeD<DIM>& in_out, const double scale) {
    for (dimension_t d = 0; d < DIM; ++d) {
        in_out[d] *= scale;
    }
}

template <dimension_t DIM>
double length(const NodeD<DIM>& node) {
    double sum = 0;
    for (dimension_t d = 0; d < DIM; ++d) {
        sum += node[d] * node[d];
    }
    return sqrt(sum);
}

/*
 * A WEB consists of "nodes" that are connected with straight "edges".
 *
 * The web is created as follows:
 * - Creates a series of edges with equal length that are connected
 * - Each edge consists of small boxes (slightly overlapping).
 */
template <dimension_t DIM>
auto CreateDataWEB = [](auto& points,
                        size_t num_entities,
                        std::uint32_t seed,
                        double world_minimum,
                        double world_maximum,
                        double box_length,
                        auto set_coordinate) {
    std::default_random_engine random_engine{seed};
    std::uniform_real_distribution<> distribution(world_minimum, world_maximum);
    std::uniform_real_distribution<> dir_dist(-1, 1);
    const double world_length = world_maximum - world_minimum;
    const int EDGE_LENGTH = world_length / 4;
    const int BOXES_PER_EDGE = EDGE_LENGTH / box_length;

    // Creates nodes
    using Node = std::array<double, DIM>;

    size_t id = 0;
    Node start;
    for (dimension_t d = 0; d < DIM; ++d) {
        start[d] = distribution(random_engine);
    }
    while (id < num_entities) {
        Node direction{};
        for (dimension_t d = 0; d < DIM; ++d) {
            direction[d] = distribution(random_engine);
        }

        double s = EDGE_LENGTH / length(direction);
        scale(direction, s);

        // constrain to world
        for (dimension_t d = 0; d < DIM; ++d) {
            double x = start[d] + direction[d];
            if (x < world_minimum || x > world_maximum) {
                direction[d] = -direction[d];
            }
        }

        // create end point
        Node end{start};
        add(end, direction);

        // Create edge consisting of boxes
        Node delta = direction;  // copy!
        scale(delta, 1. / BOXES_PER_EDGE);

        Node current = start;  // copy!
        for (int i = 0; i < BOXES_PER_EDGE && id < num_entities; ++i) {
            auto& p = points[id++];
            add(current, delta);  // TODO doing this incrementally is numerically questionable
            for (dimension_t d = 0; d < DIM; ++d) {
                double x = current[d];
                set_coordinate(p, d, x);
            }
            //std::cout << "Box: " << p << " " << start << "/" << end << " d=" << delta << std::endl;
        }

        // move on
        start = end;
    }
};

auto CreateDuplicates =
    [](auto& points, int num_unique_entries, size_t num_total_entities, std::uint32_t seed) {
        std::default_random_engine random_engine{seed};
        std::uniform_int_distribution<> distribution(0, num_unique_entries);
        for (size_t i = num_unique_entries; i < num_total_entities; ++i) {
            // copy some random other point or box
            points[i] = points[distribution(random_engine)];
        }
    };
}  // namespace

enum TestGenerator { WEB = 3, CUBE = 4, CLUSTER = 7 };

template <dimension_t DIM>
auto CreatePointDataMinMax = [](auto& points,
                                TestGenerator test_generator,
                                size_t num_entities,
                                int seed,
                                double world_minimum,
                                double world_maximum,
                                double fraction_of_duplicates) {
    auto set_coordinate_lambda = [](auto& p, dimension_t dim, auto value) {
        p[dim] = static_cast<typename std::remove_reference_t<decltype(p[0])>>(value);
    };
    // Create at least 1 unique point
    // Note that the following point generator is likely, but not guaranteed, to created unique
    // points.
    int num_unique_entries =
        static_cast<int>(1 + (num_entities - 1) * (1. - fraction_of_duplicates));
    points.reserve(num_entities);
    switch (test_generator) {
    case CUBE:
        CreateDataCUBE<DIM>(
            points, num_unique_entries, seed, world_minimum, world_maximum, set_coordinate_lambda);
        break;
    case CLUSTER:
        CreateDataCLUSTER<DIM>(
            points, num_unique_entries, seed, world_minimum, world_maximum, set_coordinate_lambda);
        break;
    default:
        assert(false);
    }

    // Create duplicates
    CreateDuplicates(points, num_unique_entries, num_entities, seed);
};

template <dimension_t DIM>
auto CreateBoxDataMinMax = [](auto& points,
                              TestGenerator test_generator,
                              size_t num_entities,
                              int seed,
                              double world_minimum,
                              double world_maximum,
                              double box_length,
                              double fraction_of_duplicates) {
    auto set_coordinate_lambda = [box_length](auto& p, dimension_t dim, auto value) {
        p.min()[dim] = value;
        p.max()[dim] = value + box_length;
    };
    // Create at least 1 unique point
    // Note that the following point generator is likely, but not guaranteed, to created unique
    // points.
    int num_unique_entries =
        static_cast<int>(1 + (num_entities - 1) * (1. - fraction_of_duplicates));
    points.reserve(num_entities);
    switch (test_generator) {
    case CUBE:
        CreateDataCUBE<DIM>(
            points, num_unique_entries, seed, world_minimum, world_maximum, set_coordinate_lambda);
        break;
    case CLUSTER:
        CreateDataCLUSTER<DIM>(
            points, num_unique_entries, seed, world_minimum, world_maximum, set_coordinate_lambda);
        break;
    case WEB:
        CreateDataWEB<DIM>(
            points,
            num_unique_entries,
            seed,
            world_minimum,
            world_maximum,
            box_length,
            set_coordinate_lambda);
        break;
    default:
        assert(false);
    }

    // Create duplicates
    CreateDuplicates(points, num_unique_entries, num_entities, seed);
};

template <dimension_t DIM>
auto CreatePointData = [](auto& points,
                          TestGenerator test_generator,
                          size_t num_entities,
                          int seed,
                          double world_length,
                          double fraction_of_duplicates = 0.) {
    CreatePointDataMinMax<DIM>(
        points, test_generator, num_entities, seed, 0, world_length, fraction_of_duplicates);
};

template <dimension_t DIM>
auto CreateBoxData = [](auto& points,
                        TestGenerator test_generator,
                        size_t num_entities,
                        int seed,
                        double world_length,
                        double box_length,
                        double fraction_of_duplicates = 0.) {
    CreateBoxDataMinMax<DIM>(
        points,
        test_generator,
        num_entities,
        seed,
        0,
        world_length,
        box_length,
        fraction_of_duplicates);
};

}  // namespace improbable::phtree::phbenchmark

#endif  // PHTREE_BENCHMARK_UTIL_H
