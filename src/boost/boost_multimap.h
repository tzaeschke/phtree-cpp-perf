// SPDX-FileCopyrightText: 2022 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef BOOST_MULTIMAP_H
#define BOOST_MULTIMAP_H

#include <boost/geometry/algorithms/comparable_distance.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/strategies/strategies.hpp>
// #include <boost/geometry/geometry.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/geometries/register/box.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/predicates.hpp>

#include <boost/geometry/geometries/adapted/std_array.hpp>

#include "boost/geometry/index/rtree.hpp"
#include "phtree/common/common.h"
#include "phtree/converter.h"
#include "phtree/filter.h"
#include <iostream>
#include <unordered_set>

#define USE_STD_ARRAY
#define USE_PH_BOX

#ifdef USE_STD_ARRAY
BOOST_GEOMETRY_REGISTER_STD_ARRAY_CS(cs::cartesian)
#ifdef USE_PH_BOX
namespace ph = improbable::phtree;
BOOST_GEOMETRY_REGISTER_BOX(ph::PhBoxD<3>, ph::PhPointD<3>, min(), max())
#endif  // USE_PH_BOX
#endif

// namespace boost { namespace geometry { namespace traits
//{
//
// template <typename T, std::size_t D> struct tag<std::array<T, D>> { using type = point_tag; };
// template <typename T, std::size_t D> struct dimension<std::array<T, D>> : boost::mpl::int_<D> {};
// template <typename T, std::size_t D> struct coordinate_type<std::array<T, D>> { using type = T;
// }; template <typename T, std::size_t D> struct coordinate_system<std::array<T, D>> { using type =
// boost::geometry::cs::cartesian; };
//
// template <typename T, std::size_t D, std::size_t Index>
// struct access<std::array<T, D>, Index> {
//     static_assert(Index < D, "Out of range");
//     using Point = std::array<T, D>;
//     using CoordinateType = typename coordinate_type<Point>::type;
//     static inline CoordinateType get(Point const& p) { return p[Index]; }
//     static inline void set(Point& p, CoordinateType const& value) { p[Index] = value; }
// };
//
// }}} // namespace boost::geometry::traits

namespace b {

using namespace boost::geometry::index;
namespace pht = improbable::phtree;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
#ifdef USE_STD_ARRAY
using point_car = std::array<double, 3>;
#else
using point_car = bg::model::point<double, 3, bg::cs::cartesian>;
#endif
#ifdef USE_PH_BOX
namespace ph = improbable::phtree;
using box_car = ph::PhBoxD<3>;
#else
using box_car = bg::model::box<point_car>;
#endif  // USE_PH_BOX

namespace {

#ifdef USE_STD_ARRAY
template <pht::dimension_t DIM, typename SCALAR>
struct PointConverter {
    using Scalar = SCALAR;
    const point_car& pre(const pht::PhPoint<DIM, SCALAR>& point) const {
        return point;
    }

    const pht::PhPoint<DIM, SCALAR>& post(const point_car& point) const {
        return point;
    }

    box_car pre_query(const pht::PhBox<DIM, SCALAR>& query_box) const {
        return {pre(query_box.min()), pre(query_box.max())};
    }

    box_car pre_query(const pht::PhPoint<DIM, SCALAR>& query_point) const {
        return {pre(query_point), pre(query_point)};
    }
};
#else
template <pht::dimension_t DIM, typename SCALAR>
struct PointConverter {
    point_car pre(const pht::PhPoint<DIM, SCALAR>& key) const {
        return {key[0], key[1], key[2]};
    }

    pht::PhPoint<DIM, SCALAR> post(const point_car& p) const {
        return {p.get<0>(), p.get<1>(), p.get<2>()};
    }

    box_car pre_query(const pht::PhBox<DIM, SCALAR>& query_box) const {
        return {pre(query_box.min()), pre(query_box.max())};
    }
};
#endif

#ifdef USE_PH_BOX
template <pht::dimension_t DIM, typename SCALAR>
struct BoxConverter {
    using Scalar = SCALAR;
    using PhBox = pht::PhBox<DIM, SCALAR>;

    const box_car& pre(const PhBox& box) const {
        return box;
    }

    PhBox post(const box_car& r) const {
        return r;
    }

    box_car pre_query(const PhBox& box) const {
        return box;
    }
};
#elif defined(USE_STD_ARRAY)
template <pht::dimension_t DIM, typename SCALAR = double>
struct BoxConverter {
    using PhBox = pht::PhBox<DIM, SCALAR>;

    box_car pre(const PhBox& box) const {
        point_car lo{box.min()[0], box.min()[1], box.min()[2]};
        point_car hi{box.max()[0], box.max()[1], box.max()[2]};
        return {lo, hi};
    }

    PhBox post(const box_car& r) const {
        auto& rlo = r.min_corner();
        auto& rhi = r.max_corner();
        pht::PhBoxD<DIM> box{rlo, rhi};
        return box;
    }

    box_car pre_query(const PhBox& box) const {
        point_car lo{box.min()[0], box.min()[1], box.min()[2]};
        point_car hi{box.max()[0], box.max()[1], box.max()[2]};
        return {lo, hi};
    }
};
#else
template <pht::dimension_t DIM, typename SCALAR>
struct BoxConverter {
    using PhBox = pht::PhBox<DIM, SCALAR>;
    using QueryBox = PhBox;

    box_car pre(const PhBox& in) const {
        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(in);
        point_car lo{box.min()[0], box.min()[1], box.min()[2]};
        point_car hi{box.max()[0], box.max()[1], box.max()[2]};
        return {lo, hi};
    }

    PhBox post(const box_car& in) const {
        auto& rlo = in.min_corner();
        auto& rhi = in.max_corner();
        pht::PhPointD<DIM> lo{rlo.get<0>(), rlo.get<1>(), rlo.get<2>()};
        pht::PhPointD<DIM> hi{rhi.get<0>(), rhi.get<1>(), rhi.get<2>()};
        pht::PhBoxD<DIM> box{lo, hi};
        return box;
    }

    auto pre_query(const QueryBox& query_box) const {
        return pre(query_box);
    }
};
#endif

template <typename PHTREE>
class IteratorBase {
  public:
    friend bool operator==(
        const IteratorBase<PHTREE>& left, const IteratorBase<PHTREE>& right) noexcept {
        return left.current_value_ptr_ == right.current_value_ptr_;
    }

    friend bool operator!=(
        const IteratorBase<PHTREE>& left, const IteratorBase<PHTREE>& right) noexcept {
        return left.current_value_ptr_ != right.current_value_ptr_;
    }
};

template <typename ITER, typename PHTREE>
class IteratorNormal {
    friend PHTREE;

  public:
    template <typename ITER_PH>
    IteratorNormal(ITER_PH&& iter_ph) noexcept : iter_{std::forward<ITER_PH>(iter_ph)} {}

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
        return iter_->second;
    }

    auto* operator->() const noexcept {
        return &iter_->second;
    }

    friend bool operator==(
        const IteratorNormal<ITER, PHTREE>& left,
        const IteratorNormal<ITER, PHTREE>& right) noexcept {
        return left.iter_ == right.iter_;
    }

    friend bool operator!=(
        const IteratorNormal<ITER, PHTREE>& left,
        const IteratorNormal<ITER, PHTREE>& right) noexcept {
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
 * The main wrapper class.
 */
template <
    pht::dimension_t DIM,
    typename T,
    typename CONVERTER = PointConverter<DIM, double>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = pht::QueryPoint>
class PhTreeMultiMap {
    using Scalar = typename CONVERTER::Scalar;
    using Key =
        typename std::conditional_t<POINT_KEYS, pht::PhPoint<DIM, Scalar>, pht::PhBox<DIM, Scalar>>;
    using PHTREE = PhTreeMultiMap<DIM, T, CONVERTER, POINT_KEYS, DEFAULT_QUERY_TYPE>;

    using Geom = typename std::conditional_t<POINT_KEYS, point_car, box_car>;
    using Entry = std::pair<Geom, T>;
    // For queries the best node capacity appears to be ¨rstar" around 16 (20 is also good).
    // THis was tested with WEB data with points and boxes.
    // TODO test insert()
    // TODO test other strategies
    using TREE = rtree<Entry, bgi::rstar<16>>; // 9; 16; 20; 25;
    using ITER = decltype(TREE().qend());

    friend IteratorBase<PHTREE>;

  public:
    using QueryBox = pht::PhBox<DIM, Scalar>;
    using KeyInternal = Geom;

    explicit PhTreeMultiMap(CONVERTER converter = CONVERTER()) : tree_{}, converter_{converter} {}

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept = default;

    void emplace(const Key& key, const T& id) {
        if constexpr (DIM == 3) {
            tree_.insert(std::make_pair(converter_.pre(key), id));
        } else {
            assert(false);
        }
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
        // TODO? This is not publicly documented....??
        //        return tree_.count(converter_.pre_query(key));

        size_t n = 0;
        auto it = tree_.qbegin(bgi::covered_by(converter_.pre_query(key)));
        for (; it != tree_.qend(); ++it) {
            ++n;
        }
        return n;
    }

    //    template <typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    //    size_t estimate_count(QueryBox query_box, QUERY_TYPE query_type = QUERY_TYPE()) const {
    //        size_t n = 0;
    //        auto counter_lambda = [&](const Key&, const BUCKET& bucket) { n += bucket.size(); };
    //        tree_->for_each(query_type(converter_.pre_query(query_box)), counter_lambda);
    //        return n;
    //    }

    auto find(const Key& key) {
        auto it = tree_.template qbegin(
            bgi::covered_by(converter_.pre_query(key)) &&
            bgi::satisfies([&](auto const& v) { return v.first == key; }));
        return IteratorNormal<ITER, PHTREE>(it);
    }

    auto find(const Key& key, const T& value) {
        auto query_box = converter_.pre_query(key);
        auto it = tree_.qbegin(bgi::covered_by(query_box) && bgi::satisfies([&](auto const& v) {
                                   return v.second == value;
                               }));
        return IteratorNormal<ITER, PHTREE>(std::move(it));
    }

    size_t erase(const Key& key, const T& value) {
        return tree_.remove(std::make_pair(converter_.pre(key), value));
    }

    template <typename ITERATOR>
    size_t erase(const ITERATOR& iterator) {
        if (iterator.iter_ == tree_.qend()) {
            return 0;
        }
        return tree_.remove(*iterator.iter_);
        // return tree_.remove(iterator.iter_, iterator.iter_);
    }

    template <typename T2>
    size_t relocate(const Key& old_key, const Key& new_key, T2&& value, bool count_equals = true) {
        erase(old_key, value);
        insert(new_key, value);
        return 1;
    }

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

    //    auto relocate_all(const Key& old_key, const Key& new_key) {
    //        return tree_->relocate(old_key, new_key);
    //    }

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

    template <typename CALLBACK, typename FILTER = pht::FilterNoOp>
    void for_each(QueryBox query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        auto predicate =
            bgi::intersects(converter_.pre_query(query_box)) && bgi::satisfies([&](auto const& v) {
                // KeyInternal k{};  // TODO
                KeyInternal k = converter_.post(v.first);
                return filter.IsBucketEntryValid(k, v.second);
            });
        auto it = tree_.qbegin(predicate);
        //        auto it = tree_.qbegin(bgi::intersects(to_region(query_box)));

        for (; it != tree_.qend(); ++it) {
            callback(converter_.post(it->first), it->second);
        }
    }

    //        template <typename FILTER = pht::FilterNoOp>
    //        auto begin(FILTER&& filter = FILTER()) const {
    //            //return CreateIterator(tree_->begin(std::forward<FILTER>(filter)));
    //            return IteratorNormal<ITER, PHTREE>(tree_.begin());
    //        }

    template <typename FILTER = pht::FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        auto it = tree_.qbegin(
            bgi::intersects(converter_.pre_query(query_box)) && bgi::satisfies([&](auto const& v) {
                KeyInternal k = converter_.post(v.first);  // TODO?
                auto id = v.second;
                return filter.IsBucketEntryValid(k, id);
            }));

        return IteratorNormal<ITER, PHTREE>(it);
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
        return IteratorNormal<ITER, PHTREE>(tree_.qend());
    }

    void clear() {
        tree_.clear();
    }

    [[nodiscard]] size_t size() const {
        return tree_.size();
    }

    [[nodiscard]] bool empty() const {
        return tree_.empty();
    }

    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    TREE tree_;
    CONVERTER converter_;
};

template <pht::dimension_t DIM, typename T, typename CONVERTER = PointConverter<DIM, double>>
using PhTreeMultiMapD = PhTreeMultiMap<DIM, T, CONVERTER>;

template <pht::dimension_t DIM, typename T, typename CONVERTER = PointConverter<DIM, float>>
using PhTreeMultiMapF = PhTreeMultiMap<DIM, T, CONVERTER>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX>
using PhTreeMultiMapBox = PhTreeMultiMap<DIM, T, CONVERTER_BOX, false, pht::QueryIntersect>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX = BoxConverter<DIM, double>>
using PhTreeMultiMapBoxD = PhTreeMultiMapBox<DIM, T, CONVERTER_BOX>;

template <pht::dimension_t DIM, typename T, typename CONVERTER_BOX = BoxConverter<DIM, float>>
using PhTreeMultiMapBoxF = PhTreeMultiMapBox<DIM, T, CONVERTER_BOX>;

}  // namespace b

#endif  // BOOST_MULTIMAP_H
