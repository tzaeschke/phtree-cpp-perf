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

#ifndef MCXME_PHTREE_H
#define MCXME_PHTREE_H

#include "PHTree.h"
#include "phtree/common/common.h"
#include "phtree/converter.h"
#include "phtree/filter.h"

namespace mcxme {

using namespace improbable::phtree;

const uint DIM = 3;
// using dimension_t = std::uint32_t;
using scalar_64_t = std::int64_t;
// using QueryPoint = PhPoint<DIM, double>;

// TODO 32 doesn´t work!
const unsigned int WIDTH = 64;
const unsigned int bitLength = WIDTH;  // 64; //4;

template <dimension_t DimInternal>
class IteratorBase {
  public:
    explicit IteratorBase() noexcept : entry_{}, finished_{true} {}

    //       explicit IteratorBase(const EntryT *current_entry) noexcept:
    //       current_entry_{current_entry} {}

    inline auto& operator*() const noexcept {
        assert(!finished_);
        return entry_;  // current_entry_->GetValue();
    }

    inline auto* operator->() const noexcept {
        assert(!finished_);
        return &entry_;  // current_entry_->GetValue();
    }

    inline friend bool operator==(
        const IteratorBase<DimInternal>& left, const IteratorBase<DimInternal>& right) noexcept {
        return left.finished_ == right.finished_ && (left.finished_ || left.entry_ == right.entry_);
    }

    inline friend bool operator!=(
        const IteratorBase<DimInternal>& left, const IteratorBase<DimInternal>& right) noexcept {
        return left.finished_ != right.finished_ ||
            (!left.finished_ && left.entry_ != right.entry_);
    }

    //        auto &second() const {
    //            return current_entry_->GetValue();
    //        }

    [[nodiscard]] inline bool IsEnd() const noexcept {
        return finished_;
    }

    //        inline EntryT *GetEntry() const noexcept {
    //            return const_cast<EntryT *>(current_entry_);
    //        }

  protected:
    void SetFinished() {
        finished_ = true;
    }

    //        void SetCurrentResult(const  *current_entry) {
    //            current_entry_ = current_entry;
    //        }

    //    protected:
    //        const EntryT *current_entry_;
    Entry<DimInternal, WIDTH> entry_;
    bool finished_;
};

template <dimension_t DimInternal>
class IteratorFull : public IteratorBase<DimInternal> {
    //        static constexpr dimension_t DIM = CONVERT::DimInternal;
    //        using SCALAR = typename CONVERT::ScalarInternal;
    //        using EntryT = typename IteratorBase<T>::EntryT;

  public:
    IteratorFull(RangeQueryIterator<DimInternal, WIDTH>* it)
    : IteratorBase<DimInternal>{}, it_{it} {
        this->finished_ = false;
        FindNextElement();
    }

    IteratorFull& operator++() noexcept {
        FindNextElement();
        return *this;
    }

    IteratorFull operator++(int) noexcept {
        IteratorFull iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void FindNextElement() noexcept {
        if (it_->hasNext()) {
            this->entry_ = it_->next();
        } else {
            this->SetFinished();
        }
    }

    //        bool IsEmpty() noexcept {
    //            return !it_l.stack_size_ == 0;
    //        }

    RangeQueryIterator<DimInternal, WIDTH>* it_;
};

//    template<dimension_t DIM, unsigned int WIDTH>
//    class RangeQueryIteratorW {
//    public:
//        RangeQueryIteratorW(RangeQueryIterator<DIM, WIDTH> *it) : it_{it} {};
//
//        inline friend bool operator==(
//                const RangeQueryIteratorW<DIM, WIDTH>& left, const IteratorBase<EntryT>& right)
//                noexcept {
//            return left.current_entry_ == right.current_entry_;
//        }
//
//        inline friend bool operator!=(
//                const RangeQueryIteratorW<DIM, WIDTH>& left, const IteratorBase<EntryT>& right)
//                noexcept {
//            return left.current_entry_ != right.current_entry_;
//        }
//
//
//    private:
//        RangeQueryIterator<DIM, WIDTH> *it_;
//    };

template <dimension_t DimInternal>
using IteratorEnd = IteratorBase<DimInternal>;

/*
 * PH-Tree main class.
 * This class is a wrapper which can implement different implementations of the PH-Tree.
 * This class support only `PhPoint` coordinates with `int64_t` scalars.
 *
 * For more information please refer to the README of this project.
 */
template <
    dimension_t DIM,
    typename T = int,
    typename CONVERTER = ConverterNoOp<DIM, scalar_64_t>,
    bool IS_BOX = CONVERTER::DimInternal != CONVERTER::DimExternal>
class PhTree {
    //    friend PhTreeDebugHelper;

    using Key = typename CONVERTER::KeyExternal;
    static constexpr dimension_t DimInternal = CONVERTER::DimInternal;

    // DimInternal==DIM indicates point keys. Box keys have DimInternal==2*DIM.
    using DEFAULT_QUERY_TYPE =
        typename std::conditional<(DIM == DimInternal), QueryPoint, QueryIntersect>::type;

  public:
    // Unless specified otherwise this is just PhBox<DIM, SCALAR_EXTERNAL>
    using QueryBox = typename CONVERTER::QueryBoxExternal;

    template <typename CONV = CONVERTER>
    explicit PhTree(CONV&& converter = CONV()) : tree_{}, converter_{converter} {}

    PhTree(const PhTree& other) = delete;

    PhTree& operator=(const PhTree& other) = delete;

    PhTree(PhTree&& other) noexcept = default;

    PhTree& operator=(PhTree&& other) noexcept = default;

    ~PhTree() noexcept = default;

    /*
     *  Attempts to build and insert a key and a value into the tree.
     *
     *  @param key The key for the new entry.
     *
     *  @param args  Arguments used to generate a new value.
     *
     *  @return  A pair, whose first element points  to the possibly inserted pair,
     *           and whose second element is a bool that is true if the pair was actually inserted.
     *
     * This function attempts to build and insert a (key, value) pair into the tree. The PH-Tree is
     * effectively a map, so if an entry with the same key was already in the tree, returns that
     * entry instead of inserting a new one.
     */
    template <typename... Args>
    void emplace(const Key& key, Args&&... args) {
        ++size_;
        tree_.insert(new_entry(key, std::forward<Args>(args)...));
        // return tree_.try_emplace(converter_.pre(key), std::forward<Args>(args)...);
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
     * // Move value from key1 to key2
     * auto iter = tree.find(key1);
     * auto value = iter.second(); // The value may become invalid in erase()
     * erase(iter);
     * emplace_hint(iter, key2, value);  // the iterator can still be used as hint here
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> emplace_hint(const ITERATOR& iterator, const Key& key, Args&&... args) {
        ++size_;
        return tree_.try_emplace(iterator, converter_.pre(key), std::forward<Args>(args)...);
    }

    /*
     * See std::map::insert().
     *
     * @return a pair consisting of the inserted element (or to the element that prevented the
     * insertion) and a bool denoting whether the insertion took place.
     */
    void insert(const Key& key, const T& value) {
        ++size_;
        tree_.insert(new_entry(key, value));
        // return tree_.insert(converter_.pre(key), value);
    }

    /*
     * See emplace().
     */
    template <typename... Args>
    std::pair<T&, bool> try_emplace(const Key& key, Args&&... args) {
        ++size_;
        return tree_.try_emplace(converter_.pre(key), std::forward<Args>(args)...);
    }

    /*
     * See emplace_hint().
     */
    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> try_emplace(const ITERATOR& iterator, const Key& key, Args&&... args) {
        ++size_;
        return tree_.try_emplace(iterator, converter_.pre(key), std::forward<Args>(args)...);
    }

    /*
     * @return the value stored at position 'key'. If no such value exists, one is added to the tree
     * and returned.
     */
    T& operator[](const Key& key) {
        return tree_[converter_.pre(key)];
    }

    /*
     * Analogous to map:count().
     *
     * @return '1', if a value is associated with the provided key, otherwise '0'.
     */
    size_t count(const Key& key) const {
        return tree_.lookup(temp_entry(key, 0)).first;
        // return tree_.count(converter_.pre(key));
    }

    /*
     * Analogous to map:find().
     *
     * Get an entry associated with a k dimensional key.
     * @param key the key to look up
     * @return an iterator that points either to the associated value or to {@code end()} if the key
     * was found
     */
    auto find(const Key& key) {
        //            if constexpr (IS_BOX) {
        //                auto k2 = converter_.pre(key);
        //                auto lo = to_vector(box2.min());
        //                auto up = to_vector(box2.max());
        //                std::pair<bool, int> result = tree_.lookupHyperRect(temp_entry(key, 0));
        //                return result;
        //            } else {
        std::pair<bool, int> result = tree_.lookup(temp_entry(key, 0));
        return result;
        //            }
    }

    /*
     * See std::map::erase(). Removes any value associated with the provided key.
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    size_t erase(const Key& key) {
        return tree_.erase(converter_.pre(key));
    }

    /*
     * See std::map::erase(). Removes any entry located at the provided iterator.
     *
     * This function attempts to use the iterator to directly erase the current entry from
     * its node. However, not all iterators provide all required information so this function
     * may resort to erase|(key, value) and thus may not be faster than that.
     *
     * Currently only iterators returned by find(key) will result in faster erase.
     *
     * @return '1' if a value was found, otherwise '0'.
     */
    template <typename ITERATOR>
    size_t erase(const ITERATOR& iterator) {
        return tree_.erase(iterator);
    }

    /*
     * This function attempts to remove a 'value' from 'old_key' and reinsert it for 'new_key'.
     *
     * The function will report _success_ in the following cases:
     * - the value was removed from the old position and reinserted at the new position
     * - the position and new position refer to the same bucket.
     *
     * The function will report _failure_ in the following cases:
     * - The value was already present in the new position
     * - The value was not present in the old position
     *
     * This method will _not_ remove the value from the old position if it is already present at the
     * new position.
     *
     * @param old_key The old position
     * @param new_key The new position
     * @return '1' if the 'value' was moved, otherwise '0'.
     */
    auto relocate(const Key& old_key, const Key& new_key) {
        return tree_.relocate_if(
            converter_.pre(old_key), converter_.pre(new_key), [](const T&) { return true; });
    }

    /*
     * Relocate (move) an entry from one position to another, subject to a predicate.
     *
     * @param old_key The old position
     * @param new_key The new position
     * @param predicate The predicate is called for every value before it is relocated.
     *                  If the predicate returns 'false', the relocation is aborted.
     * @return '1' if the 'value' was moved, otherwise '0'.
     */
    template <typename PRED>
    auto relocate_if(const Key& old_key, const Key& new_key, PRED&& predicate) {
        return tree_.relocate_if(
            converter_.pre(old_key), converter_.pre(new_key), std::forward<PRED>(predicate));
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @param callback The callback function to be called for every entry that matches the query.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     */
    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        tree_.for_each(std::forward<CALLBACK>(callback), std::forward<FILTER>(filter));
    }

    /*
     * Performs a rectangular window query. The parameters are the min and max keys which
     * contain the minimum respectively the maximum keys in every dimension.
     * @param query_box The query window.
     * @param callback The callback function to be called for every entry that matches the query.
     * The callback requires the following signature: callback(const PhPointD<DIM> &, const T &)
     * @param query_type The type of query, such as QueryIntersect or QueryInclude
     * @param filter An optional filter function. The filter function allows filtering entries and
     * sub-nodes before they are returned or traversed. Any filter function must follow the
     * signature of the default 'FilterNoOp`.
     */
    template <
        typename CALLBACK,
        typename FILTER = FilterNoOp,
        typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    void for_each(
        QueryBox query_box,
        CALLBACK&& callback,
        FILTER&& filter = FILTER(),
        QUERY_TYPE query_type = QUERY_TYPE()) const {
        tree_.for_each(
            query_type(converter_.pre_query(query_box)),
            std::forward<CALLBACK>(callback),
            std::forward<FILTER>(filter));
    }

    /*
     * Iterates over all entries in the tree. The optional filter allows filtering entries and nodes
     * (=sub-trees) before returning / traversing them. By default all entries are returned. Filter
     * functions must implement the same signature as the default 'FilterNoOp'.
     *
     * @return an iterator over all (filtered) entries in the tree,
     */
    template <typename FILTER = FilterNoOp>
    auto begin(FILTER&& filter = FILTER()) const {
        return tree_.begin(std::forward<FILTER>(filter));
    }

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
    template <typename FILTER = FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box  //,
                                                // FILTER &&filter = FILTER(),
                     // QUERY_TYPE query_type = DEFAULT_QUERY_TYPE()
    ) const {
        auto box2 = converter_.pre_query(query_box);
        auto lo = to_vector(box2.min());
        auto up = to_vector(box2.max());
        // TODO include-query
        // RangeQueryIterator<DIM, WIDTH> *it = tree_.intersectionQuery(lo, up);
        return IteratorFull<DimInternal>(tree_.intersectionQuery(lo, up));
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
    template <
        typename DISTANCE,
        typename FILTER = FilterNoOp,
        // Some magic to disable this in case of box keys, i.e. if DIM != DimInternal
        dimension_t DUMMY = DIM,
        typename std::enable_if<(DUMMY == DimInternal), int>::type = 0>
    auto begin_knn_query(
        size_t min_results,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        // We use pre() instead of pre_query() here because, strictly speaking, we want to
        // find the nearest neighbors of a (fictional) key, which may as well be a box.
        return tree_.begin_knn_query(
            min_results,
            converter_.pre(center),
            std::forward<DISTANCE>(distance_function),
            std::forward<FILTER>(filter));
    }

    /*
     * @return An iterator representing the tree's 'end'.
     */
    auto end() const {
        return IteratorEnd<DimInternal>();  // tree_.end();
    }

    /*
     * Remove all entries from the tree.
     */
    void clear() {
        size_ = 0;
        tree_.clear();  // TODO
    }

    /*
     * @return the number of entries (key/value pairs) in the tree.
     */
    [[nodiscard]] size_t size() const {
        return 1000;
        // return size_;
    }

    /*
     * @return 'true' if the tree is empty, otherwise 'false'.
     */
    [[nodiscard]] bool empty() const {
        return tree_.empty();
    }

    /*
     * @return the converter associated with this tree.
     */
    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    std::vector<unsigned long> pre(const Key& key) const {
        auto k2 = converter_.pre(key);
        std::vector<unsigned long> values(k2.size());
        for (dimension_t i = 0; i < k2.size(); ++i) {
            values[i] = k2[i];
        }
        return values;
    }

    //        std::vector<unsigned long> pre_query(const Key &key) const {
    //            auto k2 = converter_.pre_query(key);
    //            std::vector<unsigned long> values(k2.length);
    //            for (dimension_t i = 0; i < k2.length; ++i) {
    //                values[i] = k2[i];
    //            }
    //            return values;
    //        }

    template <typename QP>
    std::vector<unsigned long> to_vector(const QP& qp) const {
        std::vector<unsigned long> values(qp.size());
        for (dimension_t i = 0; i < qp.size(); ++i) {
            values[i] = qp[i];
        }
        return values;
    }

    Entry<DimInternal, WIDTH> temp_entry(const Key& key, int id) const {
        return Entry<DimInternal, WIDTH>(pre(key), id);
    }

    Entry<DimInternal, WIDTH>& new_entry(const Key& key, int id) const {
        return *new Entry<DimInternal, WIDTH>(pre(key), id);
    }

    // This is used by PhTreeDebugHelper
    const auto& GetInternalTree() const {
        return tree_;
    }

    void CheckConsistencyExternal() const {
        [[maybe_unused]] size_t n = 0;
        for ([[maybe_unused]] const auto& entry : tree_) {
            ++n;
        }
        assert(n == size());
    }

    // v16::PhTreeV16<DimInternal, T, CONVERTER> tree_;
    PHTree<DimInternal, bitLength> tree_;  // = new PHTree<1, bitLength>();
    CONVERTER converter_;
    size_t size_ = 0;
};

/*
 * Floating-point `double` version of the PH-Tree.
 * This version of the tree accepts multi-dimensional keys with floating point (`double`)
 * coordinates.
 *
 * The default implementation uses a direct lossless (in terms of numeric precision) mapping from
 * 64bit double to 64bit long integer. The mapping is defined in the Converter functions.
 * Other, lossy mapping have been shown to provide somewhat better performance (due to
 * better tree structure), but this default mapping has been chosen because it is lossless.
 *
 * For more information please refer to the README of this project.
 */
template <dimension_t DIM, typename T, typename CONVERTER = ConverterIEEE<DIM>>
using PhTreeD = PhTree<DIM, T, CONVERTER>;

/*
 * Floating-point `float` version of the PH-Tree.
 * This version of the tree accepts multi-dimensional keys with floating point (`float`)
 * coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER = ConverterFloatIEEE<DIM>>
using PhTreeF = PhTree<DIM, T, CONVERTER>;

template <dimension_t DIM, typename T, typename CONVERTER_BOX>
using PhTreeBox = PhTree<DIM, T, CONVERTER_BOX>;

/**
 * A PH-Tree that uses (axis aligned) boxes as keys.
 * The boxes are defined with 64bit 'double' floating point coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER_BOX = ConverterBoxIEEE<DIM>>
using PhTreeBoxD = PhTreeBox<DIM, T, CONVERTER_BOX>;

/**
 * A PH-Tree that uses (axis aligned) boxes as keys.
 * The boxes are defined with 32bit 'float' coordinates.
 *
 * See 'PhTreeD' for details.
 */
template <dimension_t DIM, typename T, typename CONVERTER_BOX = ConverterBoxFloatIEEE<DIM>>
using PhTreeBoxF = PhTreeBox<DIM, T, CONVERTER_BOX>;

}  // namespace mcxme

#endif  // MCXME_PHTREE_H
