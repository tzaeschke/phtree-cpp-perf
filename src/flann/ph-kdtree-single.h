// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef FLANN_KD_TREE_SINGLE_H
#define FLANN_KD_TREE_SINGLE_H

#include "flann/algorithms/dist.h"
#include "flann/algorithms/kdtree_single_index.h"
#include "flann/defines.h"
// #include <flann/flann.h>
#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include <unordered_set>

namespace flann {

using scalar_64_t = std::int64_t;

namespace pht = improbable::phtree;

/*
 * The main wrapper class
 */
template <pht::dimension_t DIM, typename T = size_t>
class KDTreeSingle {
    static_assert(std::is_same_v<size_t, T>);

    // using KeyInternal = typename CONVERTER::KeyInternal;
    using SCALAR = double;
    using Key = typename pht::PhPoint<DIM, SCALAR>;

  public:
    // using KeyInternal = typename std::conditional_t<POINT_KEYS, Point, Region>;
    using KeyInternal = Key;
    using QueryBox = typename pht::PhBox<DIM, SCALAR>;

    explicit KDTreeSingle() {
        tree_ = create_tree();
    }

    KDTreeSingle(const KDTreeSingle& other) = delete;
    KDTreeSingle& operator=(const KDTreeSingle& other) = delete;
    KDTreeSingle(KDTreeSingle&& other) noexcept = default;
    KDTreeSingle& operator=(KDTreeSingle&& other) noexcept = default;
    ~KDTreeSingle() noexcept {
        delete tree_;
    };

    void emplace(const Key& key, const T& id) {
        assert(id == tree_->size());
        Matrix<SCALAR> dataset(const_cast<SCALAR*>(&key[0]), 1, DIM);
        tree_->addPoints(dataset);
    }

    void load(const std::vector<Key>& keys) {
        Matrix<SCALAR> dataset(const_cast<SCALAR*>(&keys[0][0]), keys.size(), DIM);
        // tree_->addPoints(dataset);
        tree_->buildIndex(dataset);
        is_built_ = true;
    }

    template <typename ITERATOR>
    std::pair<T&, bool> emplace_hint(const ITERATOR&, const Key& key, const T& id) {
        emplace(key, id);
    }

    void insert(const Key& key, const T& value) {
        emplace(key, value);
    }

    void try_emplace(const Key& key, const T& value) {
        emplace(key, value);
    }

    template <typename ITERATOR>
    void try_emplace(const ITERATOR& iterator, const Key& key, const T& value) {
        emplace_hint(iterator, key, value);
    }

    const size_t X = 0;  // TODO remove
    const double R_0 = 0.1;

    size_t count(const Key& key) const {
        build();
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        tree_->radiusSearch(queries, indexes, distances, R_0, params());
        return indexes[0].size();
    }

    auto find(const Key& key) {
        build();
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        tree_->radiusSearch(queries, indexes, distances, R_0, params());
        result_.clear();
        for (auto& ind : indexes) {
            result_.insert(result_.end(), ind.begin(), ind.end());
        }
        // TODO result_ = indexes;
        return result_.begin();
    }

    auto find(const Key& key, const T& value) {
        build();
        result_.clear();
        const Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        tree_->radiusSearch(queries, indexes, distances, R_0, params());
        for (auto r : indexes) {
            for (auto r2 : r) {
                if (r2 - X == value) {
                    result_.emplace_back(r2 - X);
                    return result_.begin();
                }
            }
        }
        return result_.end();
    }

    size_t erase(const Key&, const T& value) {
        size_t result = tree_->size();
        tree_->removePoint(value);
        result -= tree_->size();
        return result;
    }

    template <typename T2>
    size_t relocate(const Key& old_key, const Key& new_key, T2&& value) {
        erase(old_key, value);
        insert(new_key, value);
        return 1;
    }

    template <typename CALLBACK, typename FILTER = pht::FilterNoOp>
    void for_each(QueryBox query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        build();
        Key& min = query_box.min();
        Key& max = query_box.max();
        Key key{};
        double diameter2 = 0;
        for (size_t d = 0; d < DIM; ++d) {
            double diff = max[d] - min[d];
            diameter2 += diff * diff;
            key[d] = (min[d] + max[d]) / 2;
        }
        double radius = std::sqrt(diameter2) / 2;
        radius = std::max(radius, R_0);
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        tree_->radiusSearch(queries, indexes, distances, radius * radius, params());
        for (auto r : indexes) {
            for (auto r2 : r) {
                auto* point = tree_->getPoint(r2);
                bool match = true;
                for (size_t d = 0; d < DIM; ++d) {
                    if (point[d] < min[d] || point[d] > max[d]) {
                        match = false;
                        break;
                    }
                }
                Key k{};  // TODO
                if (match && filter.IsEntryValid(k, r2 - X)) {
                    callback(k, r2 - X);
                }
            }
        }
    }

    template <typename FILTER = pht::FilterNoOp>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        build();
        result_.clear();
        for_each(
            query_box,
            [this](const Key&, const size_t& value) { result_.emplace_back(value); },
            std::forward<FILTER>(filter));
        return result_.begin();
    }

    template <typename DISTANCE, typename FILTER = pht::FilterNoOp>
    auto begin_knn_query(
        size_t min_results,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) {
        (void)distance_function;
        (void)filter;

        build();
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&center[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        tree_->knnSearch(queries, indexes, distances, min_results, params());

        knn_result_.clear();
        for (size_t i1 = 0; i1 < indexes.size(); ++i1) {
            for (size_t i2 = 0; i2 < indexes[i1].size(); ++i2) {
                Key k{};
                // TODO disable for performance, but required for correctness
                //                auto* x = tree_->getPoint(indexes[i1][i2]);
                //                for (size_t d = 0; d < DIM; ++d) {
                //                    k[d] = x[d];
                //                }
                knn_result_.emplace_back(
                    KNNResult{k, indexes[i1][i2] - X, std::sqrt(distances[i1][i2])});
            }
        }
        return knn_result_.begin();
    }

    auto begin() {
        build();
        result_.clear();
        Key key{};
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        tree_->radiusSearch(
            queries, indexes, distances, std::numeric_limits<SCALAR>::infinity(), params());
        for (auto r : indexes) {
            for (auto r2 : r) {
                result_.emplace_back(r2);
            }
        }
        return result_.end();
    }

    auto end() const {
        return result_.end();
    }

    auto knn_end() const {
        return knn_result_.end();
    }

    void clear() {
        delete tree_;
        tree_ = create_tree();
    }

    [[nodiscard]] size_t size() const {
        // TODO we probably should rebalance() the tree after loading. One HACK could be to
        //  do this after loading via a call to size()
        //  -->
        //  The same could be used to do a bulk-build.
        return tree_->size();
    }

    [[nodiscard]] bool empty() const {
        return tree_->size() == 0;
    }

  private:
    auto* create_tree() const {
        KDTreeSingleIndexParams params{};
        params.insert_or_assign("trees", 1);
        // auto* tree = new KDTreeIndex<L2<SCALAR>>(params);
        SCALAR data[DIM]{};
        Matrix<SCALAR> dataset(const_cast<SCALAR*>(data), 0, DIM);
        auto* tree = new KDTreeSingleIndex<L2<SCALAR>>(dataset, params);
        // tree->buildIndex(dataset);
        // tree->removePoint(0);
        return tree;
    }

    const SearchParams params() const {
        return SearchParams{FLANN_CHECKS_UNLIMITED};  // TODO unsorted?
    }

    void build() {
        if (!is_built_ && tree_->size() > 0) {
            tree_->buildIndex();
            is_built_ = true;
        }
    }

    void build() const {
        const_cast<KDTreeSingle&>(*this).build();
    }

    flann::KDTreeSingleIndex<L2<SCALAR>>* tree_;  // TODO avoid using pointer
    std::vector<T> result_{};               /// Dirty Hack!!!! TODO
    struct KNNResult {
        Key first;
        size_t second;
        SCALAR distance;
    };
    std::vector<KNNResult> knn_result_{};  /// Dirty Hack!!!! TODO
    bool is_built_ = false;
};

template <pht::dimension_t DIM, typename T>
using KDTreeSingleD = KDTreeSingle<DIM, T>;

template <pht::dimension_t DIM, typename T>
using KDTreeSingleBox = KDTreeSingle<DIM, T>;

template <pht::dimension_t DIM, typename T>
using KDTreeSingleBoxD = KDTreeSingleBox<DIM, T>;

}  // namespace flann

#endif  // FLANN_KD_TREE_SINGLE_H
