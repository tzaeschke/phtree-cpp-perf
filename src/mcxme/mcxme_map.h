// SPDX-FileCopyrightText: 2022 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef MCXME_PHTREE_H
#define MCXME_PHTREE_H

#include "PHTree.h"
#include "phtree/common/common.h"
#include "phtree/converter.h"
#include "phtree/filter.h"

namespace mcxme {

using namespace improbable::phtree;

const unsigned int WIDTH = 64;
const unsigned int bitLength = WIDTH;

template <dimension_t DimInternal>
class IteratorBase {
  public:
    explicit IteratorBase() noexcept : entry_{}, finished_{true} {}

    inline auto& operator*() const noexcept {
        assert(!finished_);
        return entry_;
    }

    inline auto* operator->() const noexcept {
        assert(!finished_);
        return &entry_;
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

    [[nodiscard]] inline bool IsEnd() const noexcept {
        return finished_;
    }

  protected:
    void SetFinished() {
        finished_ = true;
    }

    Entry<DimInternal, WIDTH> entry_;
    bool finished_;
};

template <dimension_t DimInternal>
class IteratorFull : public IteratorBase<DimInternal> {
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

    RangeQueryIterator<DimInternal, WIDTH>* it_;
};

template <dimension_t DimInternal>
using IteratorEnd = IteratorBase<DimInternal>;

/*
 * Main wrapper class.
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

    template <typename... Args>
    void emplace(const Key& key, Args&&... args) {
        ++size_;
        auto e = Entry<DimInternal, WIDTH>(pre(key), std::forward<Args>(args)...);
        tree_.insert(e);
        // tree_.insert(new_entry(key, std::forward<Args>(args)...));
    }

    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> emplace_hint(const ITERATOR& iterator, const Key& key, Args&&... args) {
        ++size_;
        return tree_.try_emplace(iterator, converter_.pre(key), std::forward<Args>(args)...);
    }

    void insert(const Key& key, const T& value) {
        ++size_;
        auto e = Entry<DimInternal, WIDTH>(pre(key), value);
        tree_.insert(e);
        // tree_.insert(new_entry(key, value));
    }

    template <typename... Args>
    std::pair<T&, bool> try_emplace(const Key& key, Args&&... args) {
        ++size_;
        return tree_.try_emplace(converter_.pre(key), std::forward<Args>(args)...);
    }

    template <typename ITERATOR, typename... Args>
    std::pair<T&, bool> try_emplace(const ITERATOR& iterator, const Key& key, Args&&... args) {
        ++size_;
        return tree_.try_emplace(iterator, converter_.pre(key), std::forward<Args>(args)...);
    }

    T& operator[](const Key& key) {
        return tree_[converter_.pre(key)];
    }

    size_t count(const Key& key) const {
        return tree_.lookup(temp_entry(key, 0)).first;
    }

    auto find(const Key& key) {
        std::pair<bool, int> result = tree_.lookup(temp_entry(key, 0));
        return result;
    }

    size_t erase(const Key& key) {
        return tree_.erase(converter_.pre(key));
    }

    template <typename ITERATOR>
    size_t erase(const ITERATOR& iterator) {
        return tree_.erase(iterator);
    }

    //    auto relocate(const Key& old_key, const Key& new_key) {
    //        return tree_.relocate_if(
    //            converter_.pre(old_key), converter_.pre(new_key), [](const T&) { return true; });
    //    }
    //
    //    template <typename PRED>
    //    auto relocate_if(const Key& old_key, const Key& new_key, PRED&& predicate) {
    //        return tree_.relocate_if(
    //            converter_.pre(old_key), converter_.pre(new_key), std::forward<PRED>(predicate));
    //    }
    //
    //    template <typename CALLBACK, typename FILTER = FilterNoOp>
    //    void for_each(
    //        CALLBACK&&, FILTER&&  // filter = FILTER()
    //    ) const {
    //        // tree_.for_each(std::forward<CALLBACK>(callback), std::forward<FILTER>(filter));
    //    }

    template <
        typename CALLBACK,
        typename FILTER = FilterNoOp,
        typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    void for_each(
        QueryBox query_box,
        CALLBACK&&,
        FILTER&&,   // filter = FILTER(),
        QUERY_TYPE  // query_type = QUERY_TYPE()
    ) const {
        for (auto it = begin_query(query_box); it != end(); ++it) {
            // TODO
            // callback();
        }
    }

    template <typename FILTER = FilterNoOp>
    auto begin(FILTER&& filter = FILTER()) const {
        return tree_.begin(std::forward<FILTER>(filter));
    }

    template <typename FILTER = FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box  //,
                                                // FILTER &&filter = FILTER(),
                                                // QUERY_TYPE query_type = DEFAULT_QUERY_TYPE()
    ) const {
        auto box2 = converter_.pre_query(query_box);
        auto lo = to_vector(box2.min());
        auto up = to_vector(box2.max());
        if constexpr (IS_BOX) {
            // TODO include-query
            // RangeQueryIterator<DIM, WIDTH> *it = tree_.intersectionQuery(lo, up);
            return IteratorFull<DimInternal>(tree_.intersectionQuery(lo, up));
        } else {
            return IteratorFull<DimInternal>(tree_.rangeQuery(lo, up));
            // return end();
        }
    }

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

    auto end() const {
        return IteratorEnd<DimInternal>();  // tree_.end();
    }

    void clear() {
        size_ = 0;
        tree_.clear();  // TODO
    }

    [[nodiscard]] size_t size() const {
        return 1000;
        // return size_;
    }

    [[nodiscard]] bool empty() const {
        return tree_.empty();
    }

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

    PHTree<DimInternal, bitLength> tree_;
    CONVERTER converter_;
    size_t size_ = 0;
};

template <dimension_t DIM, typename T, typename CONVERTER = ConverterIEEE<DIM>>
using PhTreeD = PhTree<DIM, T, CONVERTER>;

template <dimension_t DIM, typename T, typename CONVERTER = ConverterFloatIEEE<DIM>>
using PhTreeF = PhTree<DIM, T, CONVERTER>;

template <dimension_t DIM, typename T, typename CONVERTER_BOX>
using PhTreeBox = PhTree<DIM, T, CONVERTER_BOX>;

template <dimension_t DIM, typename T, typename CONVERTER_BOX = ConverterBoxIEEE<DIM>>
using PhTreeBoxD = PhTreeBox<DIM, T, CONVERTER_BOX>;

template <dimension_t DIM, typename T, typename CONVERTER_BOX = ConverterBoxFloatIEEE<DIM>>
using PhTreeBoxF = PhTreeBox<DIM, T, CONVERTER_BOX>;

}  // namespace mcxme

#endif  // MCXME_PHTREE_H
