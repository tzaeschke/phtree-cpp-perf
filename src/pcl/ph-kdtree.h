// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef PCL_KD_TREE_H
#define PCL_KD_TREE_H

#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include "pcl/kdtree/kdtree.h"
#include "pcl/kdtree/kdtree_flann.h"
#include <unordered_set>

/*
 * From: https://www2.informatik.hu-berlin.de/~sprengsz/bb-tree/
 */


namespace bb {

using scalar_64_t = std::int64_t;
const unsigned int bitLength = 64;  // 4;

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
    static_assert(std::is_same_v<uint32_t, T>);

    // using KeyInternal = typename CONVERTER::KeyInternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr pht::dimension_t DimInternal = CONVERTER::DimInternal;

  public:
    // using KeyInternal = typename std::conditional_t<POINT_KEYS, Point, Region>;
    using KeyInternal = pcl::PointXYZ;
    using QueryBox = typename CONVERTER::QueryBoxExternal;

    explicit PhTreeMultiMap(CONVERTER converter = CONVERTER()) :
        cloud_(new pcl::PointCloud<pcl::PointXYZ>()), converter_{converter}, size_{0} {
        tree_ = create_tree();
    }

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept = default;

    void emplace(const Key& key, const T& id) {
        tree_->setInputCloud();
        // TODO verify: Is this a multimap?!?!?
        ++size_;  // TODO this is bad..!
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
        for (auto r: results) {
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
        if (!is_loaded_) {
            is_loaded_ = true;
            tree_->setInputCloud(cloud_);
        }
        auto results = tree_->SearchRange(to_shape(query_box.min()), to_shape(query_box.max()));
        for (auto& r: results) {
            Key k{}; // TODO
            if (filter.IsBucketEntryValid(k, r)) {
                callback(k, r);
            }
        }
    }

    template <typename FILTER = pht::FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        auto results = tree_->SearchRange(to_shape(query_box.min()), to_shape(query_box.max()));
        result_.clear();
        result_ = results;
        //std::copy_n(shape.begin(), DIM, key.begin());
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
        tree_->RebuildDelimiters();
        return size_;
    }

    [[nodiscard]] bool empty() const {
        return tree_->getCount() == 0;
    }

    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    pcl::KdTree<KeyInternal>* create_tree() const {
        auto* tree = new pcl::KdTree<KeyInternal>(3);
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

    template <bool DUMMY = POINT_KEYS>
    typename std::enable_if<DUMMY, KeyInternal>::type to_shape(const Key& key) const {
        pcl::PointXYZ point;

        point.x = key[0];
        point.y = key[1];
        point.z = key[2];

        cloud_->points.push_back(point);
        KeyInternal p = std::vector<float>(key.begin(), key.end());
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

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 == DimInternal, Key>::type from_shape(const KeyInternal& shape) const {
        // Point** p = static_cast<Point**>(shape);
        Key key;
        std::copy_n(shape.begin(), DIM, key.begin());
//        for (pht::dimension_t d = 0; d < DIM; ++d) {
//            key[d] = p.m_pCoords[d];
//        }
        return key;
    }

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

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_;
    bool is_loaded_ = false;
    pcl::KdTree<KeyInternal>* tree_;  // TODO avoid using pointer
    CONVERTER converter_;
    size_t size_;
    std::vector<T> result_{};  /// Dirty Hack!!!! TODO
};

template <pht::dimension_t DIM, typename T, typename CONVERTER = pht::ConverterIEEE<DIM>>
using PhTreeMultiMapD = PhTreeMultiMap<DIM, T, CONVERTER>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX>
using PhTreeMultiMapBox = PhTreeMultiMap<DIM, T, CONVERTER_BOX, false, pht::QueryIntersect>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX = pht::ConverterBoxIEEE<DIM>>
using PhTreeMultiMapBoxD = PhTreeMultiMapBox<DIM, T, CONVERTER_BOX>;

}  // namespace si

#endif  // PCL_KD_TREE_H
