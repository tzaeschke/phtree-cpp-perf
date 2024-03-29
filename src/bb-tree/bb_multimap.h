// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef BB_TREE_MULTIMAP_H
#define BB_TREE_MULTIMAP_H

#include "BBTree.h"
#include "phtree/common/common.h"
#include "phtree/converter.h"
#include "phtree/filter.h"
#include <unordered_set>

/*
 * Wrapper for the BB-tree, see: https://www2.informatik.hu-berlin.de/~sprengsz/bb-tree/
 */

namespace bb {

namespace pht = improbable::phtree;

/*
 * The main wrapper class
 */
template <pht::dimension_t DIM, typename T = std::uint32_t>
class PhTreeMultiMap {
    static_assert(std::is_same_v<uint32_t, T>);
    using Key = pht::PhPoint<DIM, float>;
    using KeyD = pht::PhPoint<DIM, double>;

  public:
    using KeyInternal = std::vector<float>;
    using QueryBox = pht::PhBox<DIM, float>;
    using QueryBoxD = pht::PhBox<DIM, double>;

    explicit PhTreeMultiMap() : size_{0} {
        tree_ = create_tree();
    }

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept  //= default;
    {
        delete tree_;
    }

    void emplace(const KeyD& key, const T& id) {
        tree_->InsertObject(to_shape(key), id);
        // TODO verify: Is this a multimap?!?!?
        ++size_;  // TODO this is bad..!
    }

    void emplace(const Key& key, const T& id) {
        tree_->InsertObject(to_shape(key), id);
        // TODO verify: Is this a multimap?!?!?
        ++size_;  // TODO this is bad..!
    }

    void load(const std::vector<std::vector<float>>& keys, const std::vector<std::uint32_t>& data) {
        tree_->BulkInsert(keys, data);
        // TODO verify: Is this a multimap?!?!?
        size_ += keys.size() / DIM;  // TODO this is bad..!
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

    size_t count(const Key& key) const {
        auto results = tree_->SearchRange(to_shape(key), to_shape(key));
        return results.size();
    }

    auto find(const Key& key) {
        auto results = tree_->SearchRange(to_shape(key), to_shape(key));
        result_.clear();
        result_ = results;
        return result_.begin();
    }

    auto find(const Key& key, const T& value) {
        result_.clear();
        auto results = tree_->SearchRange(to_shape(key), to_shape(key));
        for (auto r : results) {
            if (r == value) {
                result_.emplace_back(r);
                return result_.begin();
            }
        }
        return result_.end();
    }

    size_t erase(const Key& key, const T& value) {
        tree_->DeleteObject(to_shape(key));
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
        auto results = tree_->SearchRange(to_shape(query_box.min()), to_shape(query_box.max()));
        for (auto& r : results) {
            Key k{};  // TODO
            if (filter.IsBucketEntryValid(k, r)) {
                callback(k, r);
            }
        }
    }

    template <typename CALLBACK, typename FILTER = pht::FilterNoOp>
    void for_each(QueryBoxD query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        auto results = tree_->SearchRange(to_shape(query_box.min()), to_shape(query_box.max()));
        for (auto& r : results) {
            KeyD k{};  // TODO
            if (filter.IsBucketEntryValid(k, r)) {
                callback(k, r);
            }
        }
    }

    template <typename FILTER = pht::FilterNoOp>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        auto results = tree_->SearchRange(to_shape(query_box.min()), to_shape(query_box.max()));
        result_.clear();
        for (auto& r : results) {
            KeyInternal k{};  // TODO
            if (filter.IsEntryValid(k, r)) {
                result_.emplace_back(r);
            }
        }
        // result_ = results;
        // std::copy_n(shape.begin(), DIM, key.begin());
        return result_.begin();
    }

    //    template <
    //        typename DISTANCE,
    //        typename FILTER = FilterNoOp,
    //        // Some magic to disable this in case of box keys
    //        bool DUMMY = POINT_KEYS,
    //        typename std::enable_if<DUMMY, int>::type = 0>
    //    auto begin_knn_query(
    //        size_t min_results,
    //        const Key& center,
    //        DISTANCE&& distance_function = DISTANCE(),
    //        FILTER&& filter = FILTER()) const {
    //        // We use pre() instead of pre_query() here because, strictly speaking, we want to
    //        // find the nearest neighbors of a (fictional) key, which may as well be a box.
    //        return CreateIteratorKnn(tree_->begin_knn_query(
    //            min_results,
    //            converter_.pre(center),
    //            std::forward<DISTANCE>(distance_function),
    //            std::forward<FILTER>(filter)));
    //    }

    auto end() const {
        return result_.end();
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
        // tree_->RebuildDelimiters();
        return size_;
    }

    [[nodiscard]] bool empty() const {
        return tree_->getCount() == 0;
    }

  private:
    BBTree* create_tree() const {
        auto* tree = new BBTree(3);
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

    KeyInternal to_shape(const Key& key) const {
        KeyInternal p = std::vector<float>(key.begin(), key.end());
        return p;
    }

    KeyInternal to_shape(const KeyD& key) const {
        //        KeyInternal p = std::vector<float>(key.begin(), key.end());
        //        return p;
        KeyInternal p(key.size());
        for (size_t d = 0; d < key.size(); ++d) {
            p[d] = key[d];
        }
        return p;
    }

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

    Key from_shape(const KeyInternal& shape) const {
        // Point** p = static_cast<Point**>(shape);
        Key key;
        std::copy_n(shape.begin(), DIM, key.begin());
        //        for (pht::dimension_t d = 0; d < DIM; ++d) {
        //            key[d] = p.m_pCoords[d];
        //        }
        return key;
    }

    BBTree* tree_;
    size_t size_;
    std::vector<T> result_{};  /// Dirty Hack!!!! TODO
};

template <pht::dimension_t DIM, typename T>
using PhTreeMultiMapF = PhTreeMultiMap<DIM, T>;

template <pht::dimension_t DIM, typename T>
using PhTreeMultiMapBox = PhTreeMultiMap<DIM, T>;

template <pht::dimension_t DIM, typename T>
using PhTreeMultiMapBoxF = PhTreeMultiMapBox<DIM, T>;

}  // namespace bb

#endif  // BB_TREE_MULTIMAP_H
