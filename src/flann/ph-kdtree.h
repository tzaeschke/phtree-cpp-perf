// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef FLANN_KD_TREE_H
#define FLANN_KD_TREE_H

#include "flann/algorithms/dist.h"
#include "flann/algorithms/kdtree_index.h"
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
template <
    pht::dimension_t DIM,
    typename T = std::uint32_t,
    typename CONVERTER = pht::ConverterNoOp<DIM, scalar_64_t>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = pht::QueryPoint>
class PhTreeMultiMap {
    static_assert(std::is_same_v<size_t, T>);

    // using KeyInternal = typename CONVERTER::KeyInternal;
    using SCALAR = double;
    using Key = typename pht::PhPoint<DIM, SCALAR>;
    static constexpr pht::dimension_t DimInternal = CONVERTER::DimInternal;

  public:
    // using KeyInternal = typename std::conditional_t<POINT_KEYS, Point, Region>;
    using KeyInternal = Key;
    using QueryBox = typename pht::PhBox<DIM, SCALAR>;

    explicit PhTreeMultiMap() : size_{0} {
        tree_ = create_tree();
    }

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept = default;

    void emplace(const Key& key, const T& id) {
        assert(id == tree_->size());
        Matrix<SCALAR> dataset(const_cast<SCALAR*>(&key[0]), 1, DIM);
        tree_->addPoints(dataset);
    }

    template <typename ITERATOR, typename... Args>
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

    const size_t X = 1;
    const double R_0 = 0.1;

    size_t count(const Key& key) const {
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        SearchParams params{};  // TODO unsorted?
        tree_->radiusSearch(queries, indexes, distances, R_0, params);
        return indexes[0].size();
    }

    auto find(const Key& key) {
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        SearchParams params{};  // TODO unsorted?
        tree_->radiusSearch(queries, indexes, distances, R_0, params);
        result_.clear();
        for (auto& ind : indexes) {
            result_.insert(result_.end(), ind.begin(), ind.end());
        }
        // TODO result_ = indexes;
        return result_.begin();
    }

    auto find(const Key& key, const T& value) {
        result_.clear();
        const Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        SearchParams params{};  // TODO unsorted?
        tree_->radiusSearch(queries, indexes, distances, R_0, params);
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

    size_t erase(const Key& key, const T& value) {
        tree_->removePoint(value);
        --size_;   // TODO this is bad..!
        return 1;  // TODO
    }

    //    template <typename ITERATOR>
    //    size_t erase(const ITERATOR& iterator) {
    //        static_assert(
    //            std::is_convertible_v<ITERATOR*, IteratorBase<PHTREE>*>,
    //            "erase(iterator) requires an iterator argument. For erasing by key please use "
    //            "erase(key, value).");
    //        if (iterator != end()) {
    //            auto& bucket = const_cast<BUCKET&>(*iterator.GetIteratorOfPhTree());
    //            size_t old_size = bucket.size();
    //            bucket.erase(iterator.GetIteratorOfBucket());
    //            bool success = bucket.size() < old_size;
    //            if (bucket.empty()) {
    //                success &= tree_->erase(iterator.GetIteratorOfPhTree()) > 0;
    //            }
    //            size_ -= success;
    //            return success;
    //        }
    //        return 0;
    //    }

    template <typename T2>
    size_t relocate(const Key& old_key, const Key& new_key, T2&& value, bool count_equals = true) {
        erase(old_key, value);
        insert(new_key, value);
        return 1;
    }

    template <typename CALLBACK, typename FILTER = pht::FilterNoOp>
    void for_each(QueryBox query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
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
        SearchParams params{};  // TODO unsorted?
        tree_->radiusSearch(queries, indexes, distances, radius, params);
        for (auto r : indexes) {
            for (auto r2 : r) {
                auto* point = tree_->getPoint(r2);
                for (size_t d = 0; d < DIM; ++d) {
                    if (point[d] < min[d] || point[d] > max[d]) {
                        continue;
                    }
                }
                Key k{};  // TODO
                if (filter.IsEntryValid(k, r2 - X)) {
                    callback(k, r2 - X);
                }
            }
        }
    }

    template <typename FILTER = pht::FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        result_.clear();
        for_each(
            query_box,
            [this](const Key&, const size_t& value) { result_.emplace_back(value); },
            std::forward<FILTER>(filter));
        return result_.begin();
    }

    template <
        typename DISTANCE,
        typename FILTER = pht::FilterNoOp,
        // Some magic to disable this in case of box keys
        bool DUMMY = POINT_KEYS,
        typename std::enable_if<DUMMY, int>::type = 0>
    auto begin_knn_query(
        size_t min_results,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) {
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&center[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        SearchParams params{FLANN_CHECKS_UNLIMITED};  // TODO unsorted?
        tree_->knnSearch(queries, indexes, distances, min_results, params);

        knn_result_.clear();
        for (size_t i1 = 0; i1 < indexes.size(); ++i1) {
            for (size_t i2 = 0; i2 < indexes[i1].size(); ++i2) {
                Key k{};
                // TODO
                auto* x = tree_->getPoint(indexes[i1][i2]);
                for (size_t d = 0; d < DIM; ++d) {
                    k[d] = x[d];
                }
                knn_result_.emplace_back(
                    KNNResult{k, indexes[i1][i2] - X, std::sqrt(distances[i1][i2])});
            }
        }
        return knn_result_.begin();
    }

    auto begin() {
        result_.clear();
        Key key{};
        Matrix<SCALAR> queries(const_cast<SCALAR*>(&key[0]), 1, DIM);
        std::vector<std::vector<size_t>> indexes{};
        std::vector<std::vector<SCALAR>> distances{};
        SearchParams params{};  // TODO unsorted?
        tree_->radiusSearch(
            queries, indexes, distances, std::numeric_limits<SCALAR>::infinity(), params);
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
        size_ = 0;
    }

    [[nodiscard]] size_t size() const {
        // TODO we probably should rebalance() the tree after loading. One HACK could be to
        //  do this after loading via a call to size()
        //  -->
        //  The same could be used to do a bulk-build.
        // tree_->size();
        // assert(tree_->size() == size_);
        return tree_->size();
    }

    [[nodiscard]] bool empty() const {
        return tree_->size() == 0;
    }

    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    auto* create_tree() const {
        IndexParams params{};
        params.insert_or_assign("trees", 1);
        auto* tree = new KDTreeIndex<L2<SCALAR>>(params);
        SCALAR data[DIM]{};
        Matrix<SCALAR> dataset(const_cast<SCALAR*>(data), 1, DIM);
        tree->buildIndex(dataset);
        tree->removePoint(0);
        return tree;
    }

    //    Region query_to_region(const pht::PhBoxD<DIM>& box) const {
    //        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
    //        //std::cout << "q: " << r.
    //        return r;
    //    }
    //
    //    template <bool DUMMY = POINT_KEYS>
    //    typename std::enable_if<!DUMMY, Region>::type to_shape(const Key& key) const {
    //        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(key);
    //        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
    //        return r;
    //    }

    //    template <bool DUMMY = POINT_KEYS>
    //    typename std::enable_if<DUMMY, KeyInternal>::type to_shape(const Key& key) const {
    //        pcl::PointXYZ point;
    //
    //        point.x = key[0];
    //        point.y = key[1];
    //        point.z = key[2];
    //
    //        cloud_->points.push_back(point);
    //        KeyInternal p = std::vector<float>(key.begin(), key.end());
    //        return p;
    //    }

    //    template <bool DUMMY = POINT_KEYS>
    //    typename std::enable_if<!DUMMY, Region>::type to_region(const Key& key) const {
    //        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(key);
    //        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
    //        return r;
    //    }
    //
    //    template <bool DUMMY2 = POINT_KEYS>
    //    typename std::enable_if<DUMMY2 == true, Region>::type to_region(const Key& key2) const {
    //        pht::PhPointD<DIM> key = static_cast<pht::PhPointD<DIM>>(key2);
    //        Region r = Region(&*key.begin(), &*key.begin(), DIM);
    //        return r;
    //    }
    //
    //    Key from_array(const double* a) const {
    //        Key key;
    //        for (pht::dimension_t d = 0; d < DIM; ++d) {
    //            key[d] = a[d];
    //        }
    //        return key;
    //    }
    //
    //    Key from_point(const Point& p) const {
    //        return {from_array(p.m_pCoords)};
    //    }

    //    template <pht::dimension_t DIM2 = DIM>
    //    typename std::enable_if<DIM2 == DimInternal, Key>::type from_shape(const KeyInternal&
    //    shape) const {
    //        // Point** p = static_cast<Point**>(shape);
    //        Key key;
    //        std::copy_n(shape.begin(), DIM, key.begin());
    ////        for (pht::dimension_t d = 0; d < DIM; ++d) {
    ////            key[d] = p.m_pCoords[d];
    ////        }
    //        return key;
    //    }

    //    template <pht::dimension_t DIM2 = DIM>
    //    typename std::enable_if<DIM2 != DimInternal, pht::PhBoxD<DIM>>::type from_shape(
    //        IShape* shape) const {
    //        Region r;
    //        shape->getMBR(r);
    //        // PhPointD<DIM> lo{*r.m_pLow};
    //        // PhPointD<DIM> hi{*r.m_pHigh};
    //        // PhBoxD<DIM> box{r.m_pLow, r.m_pHigh};
    //        pht::PhBoxD<DIM> box;
    //        for (pht::dimension_t d = 0; d < DIM; ++d) {
    //            box.min()[d] = r.m_pLow[d];
    //            box.max()[d] = r.m_pHigh[d];
    //        }
    //        return box;
    //    }

    Matrix<float> dataset_;
    bool is_loaded_ = false;
    flann::KDTreeIndex<L2<SCALAR>>* tree_;  // TODO avoid using pointer
    CONVERTER converter_;
    size_t size_;
    std::vector<T> result_{};  /// Dirty Hack!!!! TODO
    struct KNNResult {
        Key first;
        size_t second;
        SCALAR distance;
    };
    std::vector<KNNResult> knn_result_{};  /// Dirty Hack!!!! TODO
};

template <pht::dimension_t DIM, typename T, typename CONVERTER = pht::ConverterIEEE<DIM>>
using PhTreeMultiMapD = PhTreeMultiMap<DIM, T, CONVERTER>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX>
using PhTreeMultiMapBox = PhTreeMultiMap<DIM, T, CONVERTER_BOX, false, pht::QueryIntersect>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX = pht::ConverterBoxIEEE<DIM>>
using PhTreeMultiMapBoxD = PhTreeMultiMapBox<DIM, T, CONVERTER_BOX>;

}  // namespace flann

#endif  // FLANN_KD_TREE_H
