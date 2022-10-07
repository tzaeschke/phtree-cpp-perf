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

#ifndef BOOST_MULTIMAP_H
#define BOOST_MULTIMAP_H

#include <boost/geometry/strategies/strategies.hpp>
#include <boost/geometry/algorithms/comparable_distance.hpp>
#include <boost/geometry/algorithms/equals.hpp>
//#include <boost/geometry/algorithms/dispatch/distance.hpp>

//#include <boost/geometry/geometry.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/predicates.hpp>

#include "boost/geometry/index/rtree.hpp"
#include "phtree/common/common.h"
#include "phtree/converter.h"
#include "phtree/filter.h"
#include <unordered_set>

namespace b {

// using dimension_t = std::uint32_t;
using scalar_64_t = std::int64_t;
const unsigned int bitLength = 64;  // 4;

/*
 * PH-Tree multi-map main class.
 *
 * The PhTreeMultiMap is a wrapper around a normal PH-Tree (single value per key). The wrapper uses
 * collections to store more than one value per key.
 * By default, this multi-map is backed by std::unordered_set<T>.
 *
 * The API follows mostly the std::unordered_multimap, exceptions are pointed out.
 *
 * Differences to PhTree
 * - This is a multi-map and hence follows the std::unordered_multimap rather than std::map
 * - erase() returns an iterator instead of a pairs {iterator, bool)
 * - similar to the normal PH-Tree, emplace() returns a reference to the value instead of an
 * iterator
 *
 * For more information please refer to the README of this project.
 */

using namespace boost::geometry::index;
namespace pht = improbable::phtree;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using point_car = bg::model::point<double, 3, bg::cs::cartesian>;
using box_car = bg::model::box<point_car>;

namespace {

/*
 * Base class for the internal PH-Tree multi-map iterators.
 *
 * This base class must be distinct from the other Iterator classes because it must be agnostic of
 * the types of the fields that hold iterators. If it knew about these types then we would need
 * to provide them for the ==/!= operators, which would then make it impossible to compare
 * the generic end() iterator with any specialized iterator.
 */
template <typename PHTREE>
class IteratorBase {
    //    friend PHTREE;
    //    using T = typename PHTREE::ValueType;

  public:
    //    explicit IteratorBase() noexcept : current_value_ptr_{nullptr} {}

    friend bool operator==(
        const IteratorBase<PHTREE>& left, const IteratorBase<PHTREE>& right) noexcept {
        return left.current_value_ptr_ == right.current_value_ptr_;
    }

    friend bool operator!=(
        const IteratorBase<PHTREE>& left, const IteratorBase<PHTREE>& right) noexcept {
        return left.current_value_ptr_ != right.current_value_ptr_;
    }

    //  protected:
    //    void SetFinished() noexcept {
    //        current_value_ptr_ = nullptr;
    //    }
    //
    //    void SetCurrentValue(const T* current_value_ptr) noexcept {
    //        current_value_ptr_ = current_value_ptr;
    //    }
    //
    //  private:
    //    const T* current_value_ptr_;
};

template <typename ITER, typename PHTREE>
class IteratorNormal {//}; : public IteratorBase<PHTREE> {
    friend PHTREE;

  public:
    //explicit IteratorNormal() noexcept : IteratorBase<PHTREE>(), iter_{} {}

    template <typename ITER_PH>
    IteratorNormal(ITER_PH&& iter_ph) noexcept
    : //IteratorBase<PHTREE>(),
        iter_{std::forward<ITER_PH>(iter_ph)} {}

    IteratorNormal& operator++() noexcept {
        ++iter_;
        return *this;
    }

    IteratorNormal operator++(int) noexcept {
        IteratorNormal iterator(*this);
        ++(*this);
        return iterator;
    }

    auto& operator*() const noexcept {
        // assert(iter);
        // return const_cast<T&>(*current_value_ptr_);
        return iter_->second;
    }

    auto* operator->() const noexcept {
        // assert(current_value_ptr_);
        // return const_cast<T*>(current_value_ptr_);
        return &iter_->second;
    }

    friend bool operator==(
        const IteratorNormal<ITER, PHTREE>& left, const IteratorNormal<ITER, PHTREE>& right) noexcept {
        return left.iter_ == right.iter_;
    }

    friend bool operator!=(
        const IteratorNormal<ITER, PHTREE>& left, const IteratorNormal<ITER, PHTREE>& right) noexcept {
        return left.iter_ != right.iter_;
    }

  private:
    ITER iter_;
};

template <typename ITERATOR_PH, typename PHTREE>
class IteratorKnn : public IteratorNormal<ITERATOR_PH, PHTREE> {
  public:
    template <typename ITER_PH, typename BucketIterType>
    IteratorKnn(ITER_PH&& iter_ph, BucketIterType&& iter_bucket) noexcept
    : IteratorNormal<ITER_PH, PHTREE>(
          std::forward<ITER_PH>(iter_ph), std::forward<BucketIterType>(iter_bucket)) {}

    [[nodiscard]] double distance() const noexcept {
        return this->GetIteratorOfPhTree().distance();
    }
};

}  // namespace

/*
 * The PhTreeMultiMap class.
 */
template <
    pht::dimension_t DIM,
    typename T,
    typename CONVERTER = pht::ConverterNoOp<DIM, scalar_64_t>,
    typename BUCKET = void,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = pht::QueryPoint,
    bool IS_BOX = CONVERTER::DimInternal != CONVERTER::DimExternal>
class PhTreeMultiMap {
    //static_assert(std::is_same_v<int64_t, T>);

    using KeyInternal = typename CONVERTER::KeyInternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr pht::dimension_t DimInternal = CONVERTER::DimInternal;
    using PHTREE = PhTreeMultiMap<DIM, T, CONVERTER, BUCKET, POINT_KEYS, DEFAULT_QUERY_TYPE>;
    using ValueType = T;
    // using BucketIterType = decltype(std::declval<BUCKET>().begin());
    //  using EndType = decltype(std::declval<v16::PhTreeV16<DimInternal, BUCKET,
    //  CONVERTER>>().end());

    using Entry = std::pair<box_car, T>;
    using TREE = rtree<Entry, bgi::rstar<4>>;
    using ITER = decltype(TREE().qend());

    friend IteratorBase<PHTREE>;

  public:
    using QueryBox = typename CONVERTER::QueryBoxExternal;

    explicit PhTreeMultiMap(CONVERTER converter = CONVERTER()) : tree_{}, converter_{converter} {}

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept = default;

    /*
     *  Attempts to build and insert a key and a value into the tree.
     *
     *  @param key The key for the new entry.
     *
     *  @param args  Arguments used to generate a new value.
     *
     *  @return  A pair, whose first element points to the possibly inserted pair,
     *           and whose second element is a bool that is true if the pair was actually inserted.
     *
     * This function attempts to build and insert a (key, value) pair into the tree. The PH-Tree is
     * effectively a multi-set, so if an entry with the same key/value was already in the tree, it
     * returns that entry instead of inserting a new one.
     */
    void emplace(const Key& key, const T& id) {
        if constexpr (DIM == 3) {
            //            to_shape(key);
            //            point_car lo{};
            //            point_car hi{};
            //            box_car box{lo, hi};
            static_assert(DIM != DimInternal);
            tree_.insert(std::make_pair(to_shape(key), id));
        } else {
            assert(false);
        }
    }

    /*
     * The emplace_hint() method uses an iterator as hint for insertion.
     * The hint is ignored if it is not useful or is equal to end().
     *
     * Iterators should normally not be used after the tree has been modified. As an exception to
     * this rule, an iterator can be used as hint if it was previously used with at most one call
     * to erase() and if no other modifications occurred.
     * The following is valid:
     *
     * // Move value from key1 to key2 (if you don't want to use relocate() ).
     * auto iter = tree.find(key1);
     * auto value = iter.second(); // The value may become invalid in erase()
     * erase(iter);
     * emplace_hint(iter, key2, value);  // the iterator can still be used as hint here
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> emplace_hint(const ITERATOR&, const Key& key, const T& id) {
        emplace(key, id);
    }

    /*
     * See std::unordered_multimap::insert().
     *
     * @return a pair consisting of the inserted value (or to the value that prevented the
     * insertion if the key/value already existed) and a bool denoting whether the insertion
     * took place.
     */
    void insert(const Key& key, const T& value) {
        emplace(key, value);
    }

    /*
     * See emplace().
     */
    void try_emplace(const Key& key, const T& value) {
        emplace(key, value);
    }

    /*
     * See emplace_hint().
     */
    template <typename ITERATOR>
    void try_emplace(const ITERATOR& iterator, const Key& key, const T& value) {
        emplace_hint(iterator, key, value);
    }

    /*
     * @return '1', if a value is associated with the provided key, otherwise '0'.
     */
    size_t count(const Key& key) const {
        return tree_.count(to_shape(key));
        //        class MyVisitor : public IVisitor {
        //          public:
        //            void visitNode(const INode& /* n */) override {}
        //            void visitData(const IData&) override {
        //                ++count;
        //                // std::cout << d.getIdentifier() << std::endl;
        //                //  the ID of this data entry is an answer to the query. I will just print
        //                it to
        //                //  stdout.
        //            }
        //            void visitData(std::vector<const IData*>& /* v */) override {}
        //            size_t count;
        //        };
        //        MyVisitor v{};
        //
        //        Region r = to_region(key);
        //        Point p{};
        //        r.getCenter(p);  // TODO we are just using the center here....
        //        // TODO for boxes we should also compare the shape!
        //
        //        tree_->pointLocationQuery(p, v);
        //        return v.count;
    }

    /*
     * Estimates the result count of a rectangular window query by counting the sizes of all buckets
     * that overlap with the query box. This estimate function should be much faster than a normal
     * query, especially in trees with many entries per bucket.
     *
     * @param query_box The query window.
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     */
    //    template <typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    //    size_t estimate_count(QueryBox query_box, QUERY_TYPE query_type = QUERY_TYPE()) const {
    //        size_t n = 0;
    //        auto counter_lambda = [&](const Key&, const BUCKET& bucket) { n += bucket.size(); };
    //        tree_->for_each(query_type(converter_.pre_query(query_box)), counter_lambda);
    //        return n;
    //    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @return an iterator that points either to the first value associated with the key or
     * to {@code end()} if no value was found
     */
    auto find(const Key& key) {
        //auto query_box = to_region(key);
        // bgi::covered_by(b)
        auto it = tree_.template qbegin(bgi::covered_by(to_region(key)));
        return IteratorNormal<ITER, PHTREE>(it);
        //        class MyVisitor : public IVisitor {
        //          public:
        //            MyVisitor(std::vector<T>& r) : result{r} {};
        //            void visitNode(const INode& /* n */) override {}
        //            void visitData(const IData& d) override {
        //                result.push_back(d.getIdentifier());
        //                // std::cout << d.getIdentifier() << std::endl;
        //                //  the ID of this data entry is an answer to the query. I will just print
        //                it to
        //                //  stdout.
        //            }
        //            void visitData(std::vector<const IData*>& /* v */) override {}
        //            std::vector<T>& result;
        //        };
        //        result_.clear();
        //        MyVisitor v{result_};
        //
        //        Region r = to_region(key);
        //        Point p{};
        //        r.getCenter(p);  // TODO we are just using the center here....
        //        // TODO for boxes we should also compare the shape!
        //
        //        tree_->pointLocationQuery(p, v);
        //        return v.result.begin();
        // return CreateIterator(tree_->find(converter_.pre(key)));
    }

    /*
     * See std::unordered_multimap::find().
     *
     * @param key the key to look up
     * @param value the value to look up
     * @return an iterator that points either to the associated value of the key/value pair
     * or to {@code end()} if the key/value pair was found
     */
    auto find(const Key& key, const T& value) {
        auto query_box = to_region(key);
        // bgi::covered_by(b) // TODO
        auto it =
            tree_.template qbegin(bgi::within(query_box) && bgi::satisfies([&](auto const& v) {
                                      return v.second == value;
                                  }));
        return IteratorNormal<ITER, PHTREE>(it);
        //        class MyVisitor : public IVisitor {
        //          public:
        //            MyVisitor(const T& v, std::vector<T>& r) : value{v}, result{r} {};
        //            void visitNode(const INode& /* n */) override {}
        //            void visitData(const IData& d) override {
        //                if (d.getIdentifier() == value) {
        //                    result.push_back(d.getIdentifier());
        //                }
        //                // std::cout << d.getIdentifier() << std::endl;
        //                //  the ID of this data entry is an answer to the query. I will just print
        //                it to
        //                //  stdout.
        //            }
        //            void visitData(std::vector<const IData*>& /* v */) override {}
        //            T value;
        //            std::vector<T>& result;
        //        };
        //        result_.clear();
        //        MyVisitor v{value, result_};
        //
        //        Region r = to_region(key);
        //        Point p{};
        //        r.getCenter(p);  // TODO we are just using the center here....
        //        // TODO for boxes we should also compare the shape!
        //
        //        tree_->pointLocationQuery(p, v);
        //        return v.result.begin();
    }

    /*
     * See std::unordered_multimap::erase(). Removes the provided key/value pair if it exists.
     *
     * @return '1' if the key/value pair was found, otherwise '0'.
     */
    size_t erase(const Key& key, const T& value) {
        return tree_.remove({to_shape(key), value});
        //        auto iter_outer = tree_->find(converter_.pre(key));
        //        if (iter_outer != tree_->end()) {
        //            auto& bucket = *iter_outer;
        //            auto result = bucket.erase(value);
        //            if (bucket.empty()) {
        //                tree_->erase(iter_outer);
        //            }
        //            size_ -= result;
        //            return result;
        //        }
        //        return 0;
        // return 1;  // TODO
    }

    /*
     * See std::map::erase(). Removes any entry located at the provided iterator.
     *
     * This function uses the iterator to directly erase the entry, so it is usually faster than
     * erase(key, value).
     *
     * @return '1' if a value was found, otherwise '0'.
     */
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

    /*
     * This function attempts to remove the 'value' from 'old_key' and reinsert it for 'new_key'.
     *
     * The relocate function will report _success_ in the following cases:
     * - the value was removed from the old position and reinserted at the new position
     * - the old position and new position are identical.
     *
     * The relocate function will report _failure_ in the following cases:
     * - The value was already present in the new position
     * - The value was not present in the old position
     *
     * In case of _failure_, this function guarantees that the tree remains unchanged
     * or is returned to its original state (i.e. before the function was called).
     *
     * @param old_key The old position
     * @param new_key The new position
     * @param value The value that needs to be relocated. The relocate() method used the value's
     *              '==' operator to identify the entry that should be moved.
     * @param count_equals This setting toggles whether a relocate() between two identical keys
     *              should be counted as 'success' and return '1'. The function may still return '0'
     *              in case the keys are not in the index.
     *              Background: the intuitively correct behavior is to return '1' for identical
     *              (exising) keys. However, avoiding this check can considerably speed up
     *              relocate() calls, especially when using a ConverterMultiply.
     *
     * @return '1' if a value was found and reinserted, otherwise '0'.
     */
    template <typename T2>
    size_t relocate(const Key& old_key, const Key& new_key, T2&& value, bool count_equals = true) {
        erase(old_key, value);
        insert(new_key, value);
        return 1;
    }

    /*
     * This function attempts to remove the 'value' from 'old_key' and reinsert it for 'new_key'.
     *
     * The relocate function will report _success_ in the following cases:
     * - the value was removed from the old position and reinserted at the new position
     * - the old position and new position are identical.
     *
     * The relocate function will report _failure_ in the following cases:
     * - The value was already present in the new position
     * - The value was not present in the old position
     *
     * In case of _failure_, this function guarantees that the tree remains unchanged
     * or is returned to its original state (i.e. before the function was called).
     *
     * @param old_key The old position
     * @param new_key The new position
     * @param predicate The predicate that is used for every value at position old_key to evaluate
     *             whether it should be relocated to new_key.
     * @param count_equals This setting toggles whether a relocate() between two identical keys
     *              should be counted as 'success' and return '1'. The function may still return '0'
     *              in case the keys are not in the index.
     *              Background: the intuitively correct behavior is to return '1' for identical
     *              (exising) keys. However, avoiding this check can considerably speed up
     *              relocate() calls, especially when using a ConverterMultiply.
     *
     * @return the number of values that were relocated.
     */
    //    template <typename PREDICATE>
    //    size_t relocate_if(
    //        const Key& old_key, const Key& new_key, PREDICATE&& predicate, bool count_equals =
    //        true) { auto pair = tree_->_find_or_create_two_mm(
    //            converter_.pre(old_key), converter_.pre(new_key), count_equals);
    //        auto& iter_old = pair.first;
    //        auto& iter_new = pair.second;
    //
    //        if (iter_old.IsEnd()) {
    //            assert(iter_new.IsEnd() || !iter_new->empty());  // Otherwise remove iter_new
    //            return 0;
    //        }
    //
    //        // Are we inserting in same node and same quadrant? Or are the keys equal?
    //        if (iter_old == iter_new) {
    //            assert(old_key == new_key);
    //            return 1;
    //        }
    //
    //        size_t n = 0;
    //        auto it = iter_old->begin();
    //        while (it != iter_old->end()) {
    //            if (predicate(*it) && iter_new->emplace(std::move(*it)).second) {
    //                it = iter_old->erase(it);
    //                ++n;
    //            } else {
    //                ++it;
    //            }
    //        }
    //
    //        if (iter_old->empty()) {
    //            [[maybe_unused]] auto found = tree_->erase(iter_old);
    //            assert(found);
    //        } else if (iter_new->empty()) {
    //            [[maybe_unused]] auto found = tree_->erase(iter_new);
    //            assert(found);
    //        }
    //        return n;
    //    }

    /*
     * Relocates all values from one coordinate to another.
     * Returns an iterator pointing to the relocated data (or end(), if the relocation failed).
     */
    //    auto relocate_all(const Key& old_key, const Key& new_key) {
    //        return tree_->relocate(old_key, new_key);
    //    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default, all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @param callback The callback function to be called for every entry that matches the filter.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are passed to the callback or traversed. Any filter function must
     * follow the signature of the default 'FilterNoOp`.
     * The default 'FilterNoOp` filter matches all entries.
     */
    //    template <typename CALLBACK, typename FILTER = FilterNoOp>
    //    void for_each(CALLBACK&& callback, FILTER&& filter = FILTER()) const {
    //        class MyVisitor : public IVisitor {
    //          public:
    //            void visitNode(const INode& /* n */) override {}
    //            void visitData(const IData& d) override {
    //                if (d.getIdentifier() == value) {
    //                    result.push_back(d.getIdentifier());
    //                }
    //                // std::cout << d.getIdentifier() << std::endl;
    //                //  the ID of this data entry is an answer to the query. I will just print it
    //                to
    //                //  stdout.
    //            }
    //            void visitData(std::vector<const IData*>& /* v */) override {}
    //            CALLBACK callback_;
    //            FILTER filter_;
    //        };
    //        MyVisitor v{std::forward<CALLBACK>(callback), std::forward<FILTER>(filter)};
    //
    //
    //        PhBox<DIM> box = static_cast<PhBox<DIM>>(key);
    //        Region r = Region(box.min(), box.max(), DIM);
    //        return r;
    //        tree_->intersectsWithQuery();
    //
    //        tree_->pointLocationQuery(key, v);
    //        return v.result.begin();
    //
    //
    //
    //        tree_->for_each(
    //            NoOpCallback{},
    //            WrapCallbackFilter<CALLBACK, FILTER>{
    //                std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), converter_});
    //    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param callback The callback function to be called for every entry that matches the query
     * and filter.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     * The default 'FilterNoOp` filter matches all entries.
     */
    template <typename CALLBACK, typename FILTER = pht::FilterNoOp>
    void for_each(QueryBox query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        //        using TREE = decltype(this);
        //        class MyVisitor : public IVisitor {
        //          public:
        //            MyVisitor(CALLBACK&& cb, FILTER&& f, const TREE tree)
        //            : callback_{std::forward<CALLBACK>(cb)}
        //            , filter_{std::forward<FILTER>(f)}
        //            , tree_{tree} {};
        //            void visitNode(const INode& /* n */) override {}
        //            void visitData(const IData& d) override {
        //                KeyInternal ki{};  // TODO
        //                if (filter_.IsBucketEntryValid(ki, d.getIdentifier())) {
        //                    // Key key{};  // TODO
        //                    IShape* shape;
        //                    d.getShape(&shape);
        //                    Key k = tree_->from_shape(shape);
        //                    callback_(k, d.getIdentifier());
        //                }
        //                // std::cout << d.getIdentifier() << std::endl;
        //                //  the ID of this data entry is an answer to the query. I will just print
        //                it to
        //                //  stdout.
        //            }
        //            void visitData(std::vector<const IData*>& /* v */) override {}
        //            CALLBACK callback_;
        //            FILTER filter_;
        //            const TREE tree_;
        //        };
        //        MyVisitor v{std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), this};

        //        PhBox<DIM> box = static_cast<PhBox<DIM>>(query_box);
        //        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
        //        tree_->intersectsWithQuery(r, v);
        // auto predicate =

//        auto predicate =
//            bgi::intersects(to_region(query_box)) && bgi::satisfies([&](auto const& v) {
//                KeyInternal ki{};  // TODO
//                return filter.IsBucketEntryValid(ki, v.second);
//            });

        //        for ( auto it = tree_.qbegin(bgi::nearest(pt, 3)) ; it != tree_.qend() ; ++it )
        //        {
        //            // do something with value
        //        }

        auto it = tree_.qbegin(
            bgi::intersects(to_region(query_box)) && bgi::satisfies([&](auto const& v) {
                KeyInternal ki{};  // TODO
                return filter.IsBucketEntryValid(ki, v.second);
            }));

        for (; it != tree_.qend(); ++it) {
            Key k = from_shape(it->first);
            callback(k, it->second);
        }
        //
        //
        //
        //        tree_->template for_each<NoOpCallback, WrapCallbackFilter<CALLBACK, FILTER>>(
        //            query_type(converter_.pre_query(query_box)),
        //            {},
        //            {std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), converter_});
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default, all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @return an iterator over all (filtered) entries in the tree,
     */
    //    template <typename FILTER = FilterNoOp>
    //    auto begin(FILTER&& filter = FILTER()) const {
    //        return CreateIterator(tree_->begin(std::forward<FILTER>(filter)));
    //    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     * @return Result iterator.
     */
    template <typename FILTER = pht::FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        //        using TREE = decltype(this);
        //        class MyVisitor : public IVisitor {
        //          public:
        //            MyVisitor(FILTER&& f, std::vector<T>& r, TREE tree)
        //            : filter_{std::forward<FILTER>(f)}, result{r}, tree_{tree} {};
        //            void visitNode(const INode& /* n */) override {}
        //            void visitData(const IData& d) override {
        //                // Key k{};  // TODO
        //                IShape* shape;
        //                d.getShape(&shape);
        //                Key k = tree_->from_shape(shape);
        //                auto id = d.getIdentifier();
        //                if (filter_.IsEntryValid(k, id)) {
        //                    result.push_back(id);
        //                }
        //                // std::cout << d.getIdentifier() << std::endl;
        //                //  the ID of this data entry is an answer to the query. I will just print
        //                it to
        //                //  stdout.
        //            }
        //            void visitData(std::vector<const IData*>& /* v */) override {}
        //            FILTER filter_;
        //            std::vector<T>& result;
        //            const TREE tree_;
        //        };
        //        result_.clear();
        //        MyVisitor v{std::forward<FILTER>(filter), result_, this};

        auto it = tree_.qbegin(
            bgi::intersects(to_region(query_box)) && bgi::satisfies([&](auto const& v) {
                Key k = from_shape(v.first);  // TODO
                auto id = v.second;
                return filter.IsBucketEntryValid(k, id);
            }));

        return IteratorNormal<ITER, PHTREE>(it);
        // return result_.begin();
    }

    /*
     * Locate nearest neighbors for a given point in space.
     *
     * NOTE: This method is not (currently) available for box keys.
     *
     * @param min_results number of entries to be returned. More entries may or may not be returned
     * when several entries have the same distance.
     * @param center center point
     * @param distance_function optional distance function, defaults to euclidean distance
     * @param filter optional filter predicate that excludes nodes/entries before their distance is
     * calculated.
     * @return Result iterator.
     */
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

    /*
     * @return An iterator representing the tree's 'end'.
     */
    auto end() const {
        return IteratorNormal<ITER, PHTREE>(tree_.qend());
    }

    /*
     * Remove all entries from the tree.
     */
    void clear() {
        tree_.clear();
    }

    /*
     * @return the number of entries (key/value pairs) in the tree.
     */
    [[nodiscard]] size_t size() const {
        return tree_.size();
    }

    /*
     * @return 'true' if the tree is empty, otherwise 'false'.
     */
    //    [[nodiscard]] bool empty() const {
    //        return empty();
    //    }

    /*
     * @return the converter associated with this tree.
     */
    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    //    ISpatialIndex* create_tree() const {
    //        IStorageManager* memory = StorageManager::createNewMemoryStorageManager();
    //
    //        //        StorageManager::IBuffer* file =
    //        //        StorageManager::createNewRandomEvictionsBuffer(*memory, 10, false);
    //        //        // applies a main memory random buffer on top of the persistent storage
    //        manager
    //        //        // (LRU buffer, etc can be created the same way).
    //
    //        uint32_t indexCapacity = 100;
    //        uint32_t leafCapacity = 100;
    //
    //        // Create a new, empty, RTree with dimensionality 2, minimum load 70%, using "file" as
    //        // the StorageManager and the RSTAR splitting policy.
    //        id_type indexIdentifier;
    //        // ISpatialIndex* tree = RTree::createNewRTree(*memory, 0.7, atoi(argv[3]),
    //        atoi(argv[3]),
    //        // DIM, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
    //        ISpatialIndex* tree = RTree::createNewRTree(
    //            *memory,
    //            0.7,
    //            indexCapacity,
    //            leafCapacity,
    //            DIM,
    //            SpatialIndex::RTree::RV_RSTAR,
    //            indexIdentifier);
    //
    //        return tree;
    //    }

    box_car to_region(const pht::PhBoxD<DIM>& box) const {
        point_car lo{box.min()[0], box.min()[1], box.min()[2]};
        point_car hi{box.max()[0], box.max()[1], box.max()[2]};
        box_car b{lo, hi};
        return b;
    }

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 != DimInternal, box_car>::type to_shape(const Key& key) const {
        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(key);
        point_car lo{box.min()[0], box.min()[1], box.min()[2]};
        point_car hi{box.max()[0], box.max()[1], box.max()[2]};
        box_car b{lo, hi};
        return b;
    }

//    template <pht::dimension_t DIM2 = DIM>
//    typename std::enable_if<DIM2 == DimInternal, point_car>::type to_shape(const Key& key) const {
//        point_car p{key[0], key[1], key[2]};
//        return p;
//    }

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 != DimInternal, box_car>::type to_region(const Key& key) const {
        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(key);
        point_car lo{box.min()[0], box.min()[1], box.min()[2]};
        point_car hi{box.max()[0], box.max()[1], box.max()[2]};
        box_car b{lo, hi};
        return b;
    }

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 == DimInternal, box_car>::type to_to_region(const Key& key) const {
        point_car p{key[0], key[1], key[2]};
        // TODO can we store points directly?!?!?!?!!?
        // TODO can we store points directly?!?!?!?!!?
        // TODO can we store points directly?!?!?!?!!?
        // TODO can we store points directly?!?!?!?!!?
        // TODO can we store points directly?!?!?!?!!?
        // TODO can we store points directly?!?!?!?!!?
        // point_car hi{key[0], box.max()[1], box.max()[2]};
        box_car b{p, p};
        return b;
    }

    Key from_array(const double* a) const {
        Key key;
        for (pht::dimension_t d = 0; d < DIM; ++d) {
            key[d] = a[d];
        }
        return key;
    }

    //    Key from_point(const point_car& p) const {
    //        return {from_array(p.m_pCoords)};
    //    }
    //
    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 == DimInternal, Key>::type from_shape(const point_car& p) const {
        Key key{p.get<0>(), p.get<1>(), p.get<2>()};
        return key;
    }

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 != DimInternal, pht::PhBoxD<DIM>>::type from_shape(
        const box_car& r) const {
        auto& rlo = r.min_corner();
        auto& rhi = r.max_corner();
        pht::PhPointD<DIM> lo{rlo.get<0>(), rlo.get<1>(), rlo.get<2>()};
        pht::PhPointD<DIM> hi{rhi.get<0>(), rhi.get<1>(), rhi.get<2>()};

        pht::PhBoxD<DIM> box{lo, hi};
        return box;
    }

    // This is used by PhTreeDebugHelper
    const auto& GetInternalTree() const {
        return tree_;
    }

    TREE tree_;
    CONVERTER converter_;
    std::vector<T> result_{};  /// Dirty Hack!!!! TODO
};

/**
 * A PH-Tree multi-map that uses (axis aligned) points as keys.
 * The points are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <
    pht::dimension_t DIM,
    typename T,
    typename CONVERTER = pht::ConverterIEEE<DIM>,
    typename BUCKET = void>
using PhTreeMultiMapD = PhTreeMultiMap<DIM, T, CONVERTER, BUCKET>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX, typename BUCKET = void>
using PhTreeMultiMapBox = PhTreeMultiMap<DIM, T, CONVERTER_BOX, BUCKET, false, pht::QueryIntersect>;

/**
 * A PH-Tree multi-map that uses (axis aligned) boxes as keys.
 * The boxes are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <
    pht::dimension_t DIM,
    typename T,
    typename CONVERTER_BOX = pht::ConverterBoxIEEE<DIM>,
    typename BUCKET = void>
using PhTreeMultiMapBoxD = PhTreeMultiMapBox<DIM, T, CONVERTER_BOX, BUCKET>;

}  // namespace b

#endif  // BOOST_MULTIMAP_H
