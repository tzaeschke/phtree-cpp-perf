// SPDX-FileCopyrightText: 2022 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef LIB_SI_MULTIMAP_H
#define LIB_SI_MULTIMAP_H

#include "include/spatialindex/SpatialIndex.h"
#include "phtree/common/common.h"
#include "phtree/converter.h"
#include "phtree/filter.h"
#include <unordered_set>

namespace si {

using scalar_64_t = std::int64_t;
const unsigned int bitLength = 64;  // 4;

using namespace SpatialIndex;
namespace pht = improbable::phtree;

/*
 * The main wrapper class
 */
template <
    pht::dimension_t DIM,
    typename T,
    typename CONVERTER = pht::ConverterNoOp<DIM, scalar_64_t>,
    bool POINT_KEYS = true,
    typename DEFAULT_QUERY_TYPE = pht::QueryPoint>
class PhTreeMultiMap {
    static_assert(std::is_same_v<int64_t, T>);

    // using KeyInternal = typename CONVERTER::KeyInternal;
    using Key = typename CONVERTER::KeyExternal;
    static constexpr pht::dimension_t DimInternal = CONVERTER::DimInternal;

  public:
    // using KeyInternal = typename std::conditional_t<POINT_KEYS, Point, Region>;
    using KeyInternal = IData;
    using QueryBox = typename CONVERTER::QueryBoxExternal;

    explicit PhTreeMultiMap(CONVERTER converter = CONVERTER()) : converter_{converter}, size_{0} {
        tree_ = create_tree();
    }

    PhTreeMultiMap(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap& operator=(const PhTreeMultiMap& other) = delete;
    PhTreeMultiMap(PhTreeMultiMap&& other) noexcept = default;
    PhTreeMultiMap& operator=(PhTreeMultiMap&& other) noexcept = default;
    ~PhTreeMultiMap() noexcept = default;

    void emplace(const Key& key, const T& id) {
        tree_->insertData(0, 0, to_shape(key), id);
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
        class MyVisitor : public IVisitor {
          public:
            void visitNode(const INode& /* n */) override {}
            void visitData(const IData&) override {
                ++count;
                // std::cout << d.getIdentifier() << std::endl;
                //  the ID of this data entry is an answer to the query. I will just print it to
                //  stdout.
            }
            void visitData(std::vector<const IData*>& /* v */) override {}
            size_t count;
        };
        MyVisitor v{};

        Region r = to_region(key);
        Point p{};
        r.getCenter(p);  // TODO we are just using the center here....
        // TODO for boxes we should also compare the shape!

        tree_->pointLocationQuery(p, v);
        return v.count;
    }

    //    template <typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    //    size_t estimate_count(QueryBox query_box, QUERY_TYPE query_type = QUERY_TYPE()) const {
    //        size_t n = 0;
    //        auto counter_lambda = [&](const Key&, const BUCKET& bucket) { n += bucket.size(); };
    //        tree_->for_each(query_type(converter_.pre_query(query_box)), counter_lambda);
    //        return n;
    //    }

    auto find(const Key& key) {
        class MyVisitor : public IVisitor {
          public:
            MyVisitor(std::vector<T>& r) : result{r} {};
            void visitNode(const INode& /* n */) override {}
            void visitData(const IData& d) override {
                result.push_back(d.getIdentifier());
                // std::cout << d.getIdentifier() << std::endl;
                //  the ID of this data entry is an answer to the query. I will just print it to
                //  stdout.
            }
            void visitData(std::vector<const IData*>& /* v */) override {}
            std::vector<T>& result;
        };
        result_.clear();
        MyVisitor v{result_};

        Region r = to_region(key);
        Point p{};
        r.getCenter(p);  // TODO we are just using the center here....
        // TODO for boxes we should also compare the shape!

        tree_->pointLocationQuery(p, v);
        return v.result.begin();
    }

    auto find(const Key& key, const T& value) {
        class MyVisitor : public IVisitor {
          public:
            MyVisitor(const T& v, std::vector<T>& r) : value{v}, result{r} {};
            void visitNode(const INode& /* n */) override {}
            void visitData(const IData& d) override {
                if (d.getIdentifier() == value) {
                    result.push_back(d.getIdentifier());
                }
                // std::cout << d.getIdentifier() << std::endl;
                //  the ID of this data entry is an answer to the query. I will just print it to
                //  stdout.
            }
            void visitData(std::vector<const IData*>& /* v */) override {}
            T value;
            std::vector<T>& result;
        };
        result_.clear();
        MyVisitor v{value, result_};

        Region r = to_region(key);
        Point p{};
        r.getCenter(p);  // TODO we are just using the center here....
        // TODO for boxes we should also compare the shape!

        tree_->pointLocationQuery(p, v);
        return v.result.begin();
    }

    size_t erase(const Key& key, const T& value) {
        tree_->deleteData(to_shape(key), value);
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
        using TREE = decltype(this);
        class MyVisitor : public IVisitor {
          public:
            MyVisitor(CALLBACK&& cb, FILTER&& f, const TREE tree)
            : callback_{std::forward<CALLBACK>(cb)}
            , filter_{std::forward<FILTER>(f)}
            , tree_{tree} {};
            void visitNode(const INode& /* n */) override {}
            void visitData(const IData& d) override {
                // KeyInternal ki{};  // TODO
                //  std::cout << d.getIdentifier() << std::endl;
                if (filter_.IsBucketEntryValid(d, d.getIdentifier())) {
                    // TODO Disabling coordinate conversion improves query/WEB performance by ~5%
                    //                    Key k{};
                    IShape* shape;
                    d.getShape(&shape);
                    Key k = tree_->from_shape(shape);
                    callback_(k, d.getIdentifier());
                }
                //  the ID of this data entry is an answer to the query. I will just print it to
                //  stdout.
            }
            void visitData(std::vector<const IData*>& /* v */) override {}
            CALLBACK callback_;
            FILTER filter_;
            const TREE tree_;
        };
        static MyVisitor v{std::forward<CALLBACK>(callback), std::forward<FILTER>(filter), this};

        //        PhBox<DIM> box = static_cast<PhBox<DIM>>(query_box);
        //        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
        //        tree_->intersectsWithQuery(r, v);
        //RTree::IntersectionQuery;
        tree_->intersectsWithQuery(query_to_region(query_box), v);
    }

    //    template <typename FILTER = FilterNoOp>
    //    auto begin(FILTER&& filter = FILTER()) const {
    //        return CreateIterator(tree_->begin(std::forward<FILTER>(filter)));
    //    }

    template <typename FILTER = pht::FilterNoOp, typename QUERY_TYPE = DEFAULT_QUERY_TYPE>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) {
        using TREE = decltype(this);
        class MyVisitor : public IVisitor {
          public:
            MyVisitor(FILTER&& f, std::vector<T>& r, TREE tree)
            : filter_{std::forward<FILTER>(f)}, result{r}, tree_{tree} {};
            void visitNode(const INode& /* n */) override {}
            void visitData(const IData& d) override {
                // Key k{};      // TODO
                //                IShape* shape;
                //                d.getShape(&shape);
                //                Key k = tree_->from_shape(shape);
                //                auto id = d.getIdentifier();
                //                if (filter_.IsBucketEntryValid(k, id)) {
                //                    result.push_back(id);
                //                }
                //                IShape* shape;
                //                d.getShape(&shape);
                //                Key k = tree_->from_shape(shape);
                auto id = d.getIdentifier();
                if (filter_.IsBucketEntryValid(d, id)) {
                    result.push_back(id);
                }
                // std::cout << d.getIdentifier() << std::endl;
                //  the ID of this data entry is an answer to the query. I will just print it to
                //  stdout.
            }
            void visitData(std::vector<const IData*>& /* v */) override {}
            FILTER filter_;
            std::vector<T>& result;
            const TREE tree_;
        };
        result_.clear();
        MyVisitor v{std::forward<FILTER>(filter), result_, this};

        tree_->intersectsWithQuery(query_to_region(query_box), v);

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
        return size_;
    }

    //    [[nodiscard]] bool empty() const {
    //        return empty();
    //    }

    [[nodiscard]] const CONVERTER& converter() const {
        return converter_;
    }

  private:
    ISpatialIndex* create_tree() const {
        IStorageManager* memory = StorageManager::createNewMemoryStorageManager();

        //        StorageManager::IBuffer* file =
        //        StorageManager::createNewRandomEvictionsBuffer(*memory, 10, false);
        //        // applies a main memory random buffer on top of the persistent storage manager
        //        // (LRU buffer, etc can be created the same way).

        uint32_t indexCapacity = 9;  // 4; 9; 16; 20; 100;
        uint32_t leafCapacity = 9;   // 4; 9; 16; 20; 100;

        // Create a new, empty, RTree with dimensionality 2, minimum load 70%, using "file" as
        // the StorageManager and the RSTAR splitting policy.
        id_type indexIdentifier;
        // ISpatialIndex* tree = RTree::createNewRTree(*memory, 0.7, atoi(argv[3]), atoi(argv[3]),
        // DIM, SpatialIndex::RTree::RV_RSTAR, indexIdentifier);
        ISpatialIndex* tree = RTree::createNewRTree(
            *memory,
            0.7,
            indexCapacity,
            leafCapacity,
            DIM,
            SpatialIndex::RTree::RV_RSTAR,
            indexIdentifier);

        return tree;
    }

    Region query_to_region(const pht::PhBoxD<DIM>& box) const {
        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
        //std::cout << "q: " << r.
        return r;
    }

    template <bool DUMMY = POINT_KEYS>
    typename std::enable_if<!DUMMY, Region>::type to_shape(const Key& key) const {
        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(key);
        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
        return r;
    }

    template <bool DUMMY = POINT_KEYS>
    typename std::enable_if<DUMMY, Point>::type to_shape(const Key& key) const {
        Point p = Point(&*key.begin(), DIM);
        return p;
    }

    template <bool DUMMY = POINT_KEYS>
    typename std::enable_if<!DUMMY, Region>::type to_region(const Key& key) const {
        pht::PhBoxD<DIM> box = static_cast<pht::PhBoxD<DIM>>(key);
        Region r = Region(&*box.min().begin(), &*box.max().begin(), DIM);
        return r;
    }

    template <bool DUMMY2 = POINT_KEYS>
    typename std::enable_if<DUMMY2 == true, Region>::type to_region(const Key& key2) const {
        pht::PhPointD<DIM> key = static_cast<pht::PhPointD<DIM>>(key2);
        Region r = Region(&*key.begin(), &*key.begin(), DIM);
        return r;
    }

    Key from_array(const double* a) const {
        Key key;
        for (pht::dimension_t d = 0; d < DIM; ++d) {
            key[d] = a[d];
        }
        return key;
    }

    Key from_point(const Point& p) const {
        return {from_array(p.m_pCoords)};
    }

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 == DimInternal, Key>::type from_shape(const IShape* shape) const {
        // Point** p = static_cast<Point**>(shape);
        Point p{};
        shape->getCenter(p);
        Key key;
        for (pht::dimension_t d = 0; d < DIM; ++d) {
            key[d] = p.m_pCoords[d];
        }
        return key;
    }

    template <pht::dimension_t DIM2 = DIM>
    typename std::enable_if<DIM2 != DimInternal, pht::PhBoxD<DIM>>::type from_shape(
        IShape* shape) const {
        Region r;
        shape->getMBR(r);
        // PhPointD<DIM> lo{*r.m_pLow};
        // PhPointD<DIM> hi{*r.m_pHigh};
        // PhBoxD<DIM> box{r.m_pLow, r.m_pHigh};
        pht::PhBoxD<DIM> box;
        for (pht::dimension_t d = 0; d < DIM; ++d) {
            box.min()[d] = r.m_pLow[d];
            box.max()[d] = r.m_pHigh[d];
        }
        return box;
    }

    SpatialIndex::ISpatialIndex* tree_;
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

#endif  // LIB_SI_MULTIMAP_H
