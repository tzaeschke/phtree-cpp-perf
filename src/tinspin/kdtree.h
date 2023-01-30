// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef TINSPIN_KD_TREE_H
#define TINSPIN_KD_TREE_H

#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include "src/util/ph-util.h"
#include <iostream>
#include <unordered_set>

namespace tinspin {

using namespace improbable::phtree;

namespace {

struct FilterNoOp {
    /*
     * @param key The key/coordinate of the entry.
     * @param value The value of the entry. For MultiMaps, this is a container of values.
     * @returns This default implementation always returns `true`.
     */
    template <typename KeyT, typename ValueT>
    constexpr bool IsEntryValid(const KeyT& /*key*/, const ValueT& /*value*/) const noexcept {
        return true;
    }
};

template <typename T>
struct FilterValue {
    template <typename KeyT>
    constexpr bool IsEntryValid(const KeyT& /*key*/, const T& value) const noexcept {
        return value_ == value;
    }
    const T& value_;
};

template <typename Key>
static bool isEnclosed(const Key& point, const Key& min, const Key& max) {
    for (size_t i = 0; i < point.size(); ++i) {
        if (point[i] < min[i] || point[i] > max[i]) {
            return false;
        }
    }
    return true;
}
}  // namespace

template <typename Key, typename T>
class Node;

struct KDStats {
    size_t nNodes = 0;
    size_t maxDepth = 0;
};

template <typename Key, typename T>
class KDEntryDist {
  public:
    KDEntryDist(){};
    KDEntryDist(Node<Key, T>* node, double dist) {
        entry_ = node;
        distance_ = dist;
    }

    void set(Node<Key, T>* node, double dist) {
        entry_ = node;
        distance_ = dist;
    }

    double dist() const noexcept {
        return distance_;
    }

    const Key& point() const noexcept {
        return entry_->point();
    }

    T& value() const noexcept {
        return entry_->value();
    }

    Node<Key, T>* _node() const noexcept {
        return entry_;
    }

  private:
    Node<Key, T>* entry_;
    double distance_;
};

/**
 * Node class for the kdtree.
 */
template <typename Key, typename T>
class Node {
  public:
    Node(const Key& p, const T& value, int dim) {
        coordinate_ = p;
        value_ = value;
        dim_ = dim;
    }

    ~Node() {
        if (left_ != nullptr) {
            delete left_;
        }
        if (right_ != nullptr) {
            delete right_;
        }
    }

    Node* getClosestNodeOrAddPoint(const Key& p, const T& value, int dims) {
        // Find best sub-node.
        // If there is no node, we create one and return null
        if (p[dim_] >= coordinate_[dim_]) {
            if (right_ != nullptr) {
                return right_;
            }
            right_ = new Node(p, value, (dim_ + 1) % dims);
            return nullptr;
        }
        if (left_ != nullptr) {
            return left_;
        }
        left_ = new Node(p, value, (dim_ + 1) % dims);
        return nullptr;
    }

    const Key& getKey() const noexcept {
        return coordinate_;
    }

    const T& getValue() const noexcept {
        return value_;
    }

    Node* getLo() const noexcept {
        return left_;
    }

    Node* getHi() const noexcept {
        return right_;
    }

    void setLeft(Node* left) noexcept {
        if (left_ != nullptr) {
            delete left_;
        }
        left_ = left;
    }

    void setRight(Node* right) noexcept {
        if (right_ != nullptr) {
            delete right_;
        }
        right_ = right;
    }

    void setKeyValue(const Key& key, const T& value) noexcept {
        coordinate_ = key;
        value_ = value;
    }

    const Key& point() const noexcept {
        return coordinate_;
    }

    const T& value() const noexcept {
        return value_;
    }

    void stats(KDStats& s, size_t depth) const {
        s.nNodes++;
        if (depth > s.maxDepth) {
            s.maxDepth = depth;
        }
        if (left_ != nullptr) {
            left_->stats(s, depth + 1);
        }
        if (right_ != nullptr) {
            right_->stats(s, depth + 1);
        }
    }

    void check_consistency(
        const Key& lower_bound, const Key& upper_bound, size_t& n, size_t d_parent) const {
        assert(dim_ == (d_parent + 1u) % coordinate_.size());
        for (size_t d = 0; d < lower_bound.size(); ++d) {
            assert(coordinate_[d] >= lower_bound[d]);
            assert(coordinate_[d] <= upper_bound[d]);  // TODO <- backport?
        }

        ++n;
        if (left_ != nullptr) {
            Key upper_bound_2{upper_bound};
            upper_bound_2[dim_] = std::min(upper_bound_2[dim_], coordinate_[dim_]);
            left_->check_consistency(lower_bound, upper_bound_2, n, dim_);
        }
        if (right_ != nullptr) {
            Key lower_bound_2{lower_bound};
            lower_bound_2[dim_] = std::max(lower_bound_2[dim_], coordinate_[dim_]);
            right_->check_consistency(lower_bound_2, upper_bound, n, dim_);
        }
    }

    bool isLeaf() const noexcept {
        return left_ == nullptr && right_ == nullptr;
    }

    size_t dim() const noexcept {
        return dim_;
    }

  private:
    Key coordinate_;
    T value_;
    Node* left_ = nullptr;
    Node* right_ = nullptr;
    size_t dim_;
};

template <typename Key, typename T>
class KDIteratorBase {
  public:
    explicit KDIteratorBase() noexcept : node_{nullptr} {}

    inline auto& operator*() const noexcept {
        assert(node_ != nullptr);
        return node_->getValue();
    }

    inline auto* operator->() const noexcept {
        assert(node_ != nullptr);
        return &node_->getValue();
    }

    inline friend bool operator==(
        const KDIteratorBase<Key, T>& left, const KDIteratorBase<Key, T>& right) noexcept {
        // TODO compare stack status left/right/key
        return left.node_ == right.node_;
    }

    inline friend bool operator!=(
        const KDIteratorBase<Key, T>& left, const KDIteratorBase<Key, T>& right) noexcept {
        return left.node_ != right.node_;
    }

    auto _node() const noexcept {
        return node_;
    }

  protected:
    bool IsEnd() const noexcept {
        return this->_node() == nullptr;
    }

    void SetFinished() noexcept {
        node_ = nullptr;
    }

    void SetCurrentResult(Node<Key, T>* node) noexcept {
        node_ = node;
    }

  protected:
    Node<Key, T>* node_ = nullptr;
};

/**
 * iterator.
 */
template <typename Key, typename T, typename FILTER>
class KDIterator : public KDIteratorBase<Key, T> {
    struct IteratorPos {
        Node<Key, T>* node_;
        int depth_;
        bool doLeft, doKey, doRight;

        void set(Node<Key, T>* node, const Key& min, const Key& max, size_t depth) {
            node_ = node;
            depth_ = depth;
            const Key& key = node_->getKey();
            auto dims = min.size();
            dimension_t dim = depth % dims;
            // TODO backport -> invariant problem, invariant can be
            //   broken in additional cases with coordinayte duplicates (multiple dimensions are
            //   equal).
            doLeft = min[dim] <= key[dim];
            doRight = max[dim] >= key[dim];  // TODO backport to Java !!!!!!!!!!!!!!!!!!!!
            doKey = doLeft || doRight || key[dim] == min[dim] || key[dim] == max[dim];
        }
    };

    // TODO use std::stack instead. (perf test!)
    class IteratorStack {
      public:
        IteratorStack(size_t size = 10) : stack(size) {}

        bool isEmpty() {
            return size == 0u;
        }

        IteratorPos prepareAndPush(Node<Key, T>* node, const Key& min, const Key& max, int depth) {
            if (size == stack.size()) {
                stack.emplace_back();
            }
            IteratorPos& ni = stack[size++];

            ni.set(node, min, max, depth);
            return ni;
        }

        IteratorPos& peek() {
            return stack[size - 1];
        }

        IteratorPos& pop() {
            return stack[--size];
        }

        void clear() {
            size = 0;
        }

      private:
        std::vector<IteratorPos> stack{};
        size_t size = 0;
    };

  public:
    // end()
    KDIterator() : KDIteratorBase<Key, T>(), stack_(0), min_{}, max_{} {}

    // begin_query()
    template <typename F = FilterNoOp>
    KDIterator(Node<Key, T>* root, const Key& min, const Key& max, F&& filter = F())
    : KDIteratorBase<Key, T>(), stack_{}, min_{min}, max_{max}, filter_(std::forward<F>(filter)) {
        if (root != nullptr) {
            stack_.prepareAndPush(root, min, max, 0);
            findNext();
        } else {
            this->SetFinished();
        }
    }

    KDIterator& operator++() noexcept {
        assert(!this->IsEnd());
        findNext();
        return *this;
    }

    KDIterator operator++(int) noexcept {
        assert(!this->IsEnd());
        KDIterator iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void findNext() {
        while (!stack_.isEmpty()) {
            IteratorPos& itPos = stack_.peek();
            Node<Key, T>* node = itPos.node_;
            if (itPos.doLeft && node->getLo() != nullptr) {
                itPos.doLeft = false;
                stack_.prepareAndPush(node->getLo(), min_, max_, itPos.depth_ + 1);
                continue;
            }
            if (itPos.doKey) {
                itPos.doKey = false;
                if (isEnclosed(node->getKey(), min_, max_) &&
                    filter_.IsEntryValid(node->getKey(), node->getValue())) {
                    this->SetCurrentResult(node);
                    return;
                }
            }
            if (itPos.doRight && node->getHi() != nullptr) {
                itPos.doRight = false;
                stack_.prepareAndPush(node->getHi(), min_, max_, itPos.depth_ + 1);
                continue;
            }
            stack_.pop();
        }
        this->SetFinished();
    }

  private:
    IteratorStack stack_;
    Key min_;
    Key max_;
    FILTER filter_;
};

template <typename Key, typename T>
class KDIteratorKnn : public KDIteratorBase<Key, T> {
    using Candidates = std::vector<KDEntryDist<Key, T>>;
    using CandidatesIter = decltype(Candidates{}.begin());

  public:
    KDIteratorKnn(Candidates&& result) noexcept
    : KDIteratorBase<Key, T>(), result_{std::move(result)}, iter_{result_.begin()} {
        if (iter_ != result_.end()) {
            this->SetCurrentResult(iter_->_node());
        } else {
            this->SetFinished();
        }
    }

    KDIteratorKnn& operator++() noexcept {
        assert(!this->IsEnd());
        findNext();
        return *this;
    }

    KDIteratorKnn operator++(int) noexcept {
        KDIteratorKnn iterator(*this);
        ++(*this);
        return iterator;
    }

    [[nodiscard]] double distance() const noexcept {
        return iter_->dist();
    }

    const Key& first() const noexcept {
        return iter_->point();
    }

  private:
    void findNext() {
        assert(iter_ != result_.end());
        ++iter_;
        if (iter_ != result_.end()) {
            this->SetCurrentResult(iter_->_node());
        } else {
            this->SetFinished();
        }
    }

    Candidates result_;
    CandidatesIter iter_;
};

/**
 * A simple KD-Tree implementation.
 *
 * By T. Zäschke
 */
template <typename Key, typename T>
class KDTree {
    using SCALAR = std::remove_reference_t<decltype(Key{}[0])>;
    static_assert(std::is_floating_point_v<SCALAR>);
    // TODO this should be infinity for floats.
    static constexpr SCALAR SCALAR_MAX = std::is_floating_point_v<SCALAR>
        ? std::numeric_limits<SCALAR>::infinity()
        : std::numeric_limits<SCALAR>::max();
    static constexpr SCALAR SCALAR_MIN = std::is_floating_point_v<SCALAR>
        ? -std::numeric_limits<SCALAR>::infinity()
        : std::numeric_limits<SCALAR>::min();

    using QueryBox = tinspin::Box<Key>;

    static const bool DEBUG = false;

    struct PointEntry {
        Key& first;
        T& second;
    };

    struct RemoveResult {
        Node<Key, T>* node = nullptr;
        Node<Key, T>* nodeParent = nullptr;
        double best;
        dimension_t dim;
    };

  public:
    using KeyInternal = Key;

    //  public:
    //    void main(char* args) {
    //        for (int i = 0; i < 10; i++) {
    //            test(i);
    //        }
    //    }
    //
    //  private:
    //    static void test(int r) {
    //        //		const Key&[] point_list = {{2,3}, {5,4}, {9,6}, {4,7}, {8,1}, {7,2}};
    //        const Key&[] point_list = new double[500000][14];
    //        Random R = new Random(r);
    //        for (const Key&p : point_list) {
    //            Arrays.setAll(p, [](int i) { return (double)R.nextInt(100); });
    //        }
    //        KDTree<double[]> tree = create(point_list[0].length);
    //        for (const Key&data : point_list) {
    //            tree.insert(data, data);
    //        }
    //        //	    System.out.println(tree.toStringTree());
    //        for (const Key&key : point_list) {
    //            if (!tree.containsExact(key)) {
    //                throw new IllegalStateException("" + Arrays.toString(key));
    //            }
    //        }
    //        for (const Key&key : point_list) {
    //            System.out.println(Arrays.toString(tree.queryExact(key)));
    //        }
    //        //	    System.out.println(tree.toStringTree());
    //
    //        for (const Key&key : point_list) {
    //            //			System.out.println(tree.toStringTree());
    //            System.out.println("kNN query: " + Arrays.toString(key));
    //            QueryIteratorKNN<PointEntryDist<double[]>> iter = tree.queryKNN(key, 1);
    //            if (!iter.hasNext()) {
    //                throw new IllegalStateException("kNN() failed: " + Arrays.toString(key));
    //            }
    //            const Key&answer = iter.next().point();
    //            if (answer != key && !Arrays.equals(answer, key)) {
    //                throw new IllegalStateException(
    //                    "Expected " + Arrays.toString(key) + " but got " +
    //                    Arrays.toString(answer));
    //            }
    //        }
    //
    //        for (const Key&key : point_list) {
    //            //			System.out.println(tree.toStringTree());
    //            System.out.println("Removing: " + Arrays.toString(key));
    //            if (!tree.containsExact(key)) {
    //                throw new IllegalStateException("containsExact() failed: " +
    //                Arrays.toString(key));
    //            }
    //            const Key&answer = tree.remove(key);
    //            if (answer != key && !Arrays.equals(answer, key)) {
    //                throw new IllegalStateException(
    //                    "Expected " + Arrays.toString(key) + " but got " +
    //                    Arrays.toString(answer));
    //            }
    //        }
    //    }

  public:
    KDTree(size_t dims = 3) : dims_{dims} {
        if (DEBUG) {
            std::cout << "Warning: DEBUG enabled" << std::endl;
        }
    }

    ~KDTree() {
        delete root_;
    }

    /**
     * Insert a key-value pair.
     * @param key the key
     * @param value the value
     */
    void emplace(const Key& key, const T& value) {
        size_++;
        if (root_ == nullptr) {
            root_ = new Node(key, value, 0);
            return;
        }
        Node<Key, T>* n = root_;
        while ((n = n->getClosestNodeOrAddPoint(key, value, dims_)) != nullptr)
            ;
    }

    void insert(const Key& key, const T& value) {
        size_++;
        if (root_ == nullptr) {
            root_ = new Node(key, value, 0);
            return;
        }
        Node<Key, T>* n = root_;
        while ((n = n->getClosestNodeOrAddPoint(key, value, dims_)) != nullptr)
            ;
    }

    /**
     * Check whether a given key exists.
     * @param key the key to check
     * @return true iff the key exists
     */
    size_t count(const Key& key) {
        size_t n = 0;
        for_each({key, key}, [&n](const Key&, const T&) { ++n; });
        return n;
    }

    /**
     * Get the value associates with the key.
     * @param key the key to look up
     * @return the value for the key or 'null' if the key was not found
     */
    auto find(const Key& key) {
        return begin_query({key, key});
    }

    auto find(const Key& key, const T& value) {
        // TODO implement special query on Key i.o. Box
        // return begin_query({key, key}, [&value](const Key&, const T& v) { return value == v; });
        return begin_query({key, key}, FilterValue<T>{value});
    }

  private:
    template <typename PRED>
    Node<Key, T>* findNodeExact(const Key& key, RemoveResult& resultDepth, const PRED& pred_fn) {
        if (root_ == nullptr) {
            return nullptr;
        }
        return invariantBroken ? findNodeExactSlow(key, root_, nullptr, resultDepth, pred_fn)
                               : findNodeExactFast(key, nullptr, resultDepth, pred_fn);
    }

    // TODO remove this, this can be used only once!  (or multiple times if we never have identical
    // coordinates)
    //    is id bad enough if any arbitrary dimension is equal?
    template <typename PRED>
    Node<Key, T>* findNodeExactFast(
        const Key& key, Node<Key, T>* parent, RemoveResult& resultDepth, const PRED& pred_fn) {
        Node<Key, T>* n = root_;
        do {
            const Key& nodeKey = n->getKey();
            double nodeX = nodeKey[n->dim()];
            double keyX = key[n->dim()];
            if (keyX == nodeX && key == nodeKey && pred_fn(n->getValue())) {
                resultDepth.dim = n->dim();
                resultDepth.nodeParent = parent;
                return n;
            }
            parent = n;
            n = (keyX >= nodeX) ? n->getHi() : n->getLo();
        } while (n != nullptr);
        return n;
    }

    template <typename PRED>
    Node<Key, T>* findNodeExactSlow(
        const Key& key,
        Node<Key, T>* n,
        Node<Key, T>* parent,
        RemoveResult& resultDepth,
        const PRED& pred_fn) {
        do {
            auto& nodeKey = n->getKey();
            double nodeX = nodeKey[n->dim()];
            double keyX = key[n->dim()];
            if (keyX == nodeX) {
                if (key == nodeKey && pred_fn(n->getValue())) {
                    resultDepth.dim = n->dim();
                    resultDepth.nodeParent = parent;
                    return n;
                }
                // Broken invariant? We need to check the 'lower' part as well...
                if (n->getLo() != nullptr) {  // TODO backport!!!!!! "parent" -> "n"
                    Node<Key, T>* n2 = findNodeExactSlow(key, n->getLo(), n, resultDepth, pred_fn);
                    if (n2 != nullptr) {
                        return n2;
                    }
                }
            }
            parent = n;
            n = (keyX >= nodeX) ? n->getHi() : n->getLo();
        } while (n != nullptr);
        return n;
    }

    /**
     * Remove a key.
     * @param key key to remove
     * @return the value associated with the key or 'null' if the key was not found
     */
  public:
    size_t erase(const Key& key, const T& value) {
        return erase_if(key, [&value](const T& v) { return v == value; });
    }

    size_t erase(const Key& key) {
        return erase_if(key, [](const T&) { return true; });
    }

    template <typename PRED>
    size_t erase_if(const Key& key, const PRED& pred_fn) {
        if (root_ == nullptr) {
            return 0;
        }

        invariantBroken = true;  // TODO backport

        // find
        RemoveResult removeResult{};
        Node<Key, T>* eToRemove = findNodeExact(key, removeResult, pred_fn);
        if (eToRemove == nullptr) {
            return 0;
        }

        // remove
        if (eToRemove == root_ && size_ == 1) {
            root_ = nullptr;
            size_ = 0;
            invariantBroken = false;
        }

        // find replacement
        //   removeResult.nodeParent = nullptr; // TODO backport!
        while (eToRemove != nullptr && !eToRemove->isLeaf()) {
            // recurse
            auto dim = removeResult.dim;
            removeResult.node = nullptr;
            if (eToRemove->getHi() != nullptr) {
                // get replacement from right
                // This is preferable, because it cannot break the invariant
                removeResult.best = SCALAR_MAX;
                removeMinLeaf(eToRemove->getHi(), eToRemove, dim, removeResult);
            } else if (eToRemove->getLo() != nullptr) {
                // get replacement from left
                removeResult.best = SCALAR_MIN;
                removeMaxLeaf(eToRemove->getLo(), eToRemove, dim, removeResult);
            }
            eToRemove->setKeyValue(removeResult.node->getKey(), removeResult.node->getValue());
            eToRemove = removeResult.node;
        }
        // leaf node
        Node<Key, T>* parent = removeResult.nodeParent;
        if (parent != nullptr) {
            if (parent->getLo() == eToRemove) {
                parent->setLeft(nullptr);
            } else if (parent->getHi() == eToRemove) {
                parent->setRight(nullptr);
            } else {
                assert(false);
            }
        }
        size_--;
        return 1;
    }

  private:
    void removeMinLeaf(
        Node<Key, T>* node, Node<Key, T>* parent, dimension_t dim, RemoveResult& result) {
        // Split in 'interesting' dimension
        if (dim == node->dim()) {
            // We strictly look for leaf nodes with left==null
            //  -> left!=null means the left child is at least as small as the current node
            if (node->getLo() != nullptr) {
                removeMinLeaf(node->getLo(), node, dim, result);
            } else if (node->getKey()[dim] <= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = node->getKey()[dim];
                result.dim = node->dim();
            }
        } else {
            // split in any other dimension.
            // First, check local key.
            double localX = node->getKey()[dim];
            if (localX <= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = localX;
                result.dim = node->dim();
            }
            if (node->getLo() != nullptr) {
                removeMinLeaf(node->getLo(), node, dim, result);
            }
            if (node->getHi() != nullptr) {
                removeMinLeaf(node->getHi(), node, dim, result);
            }
        }
    }

    void removeMaxLeaf(
        Node<Key, T>* node, Node<Key, T>* parent, dimension_t dim, RemoveResult& result) {
        // Split in 'interesting' dimension
        if (dim == node->dim()) {
            // We strictly look for leaf nodes with left==null
            if (node->getHi() != nullptr) {
                removeMaxLeaf(node->getHi(), node, dim, result);
            } else if (node->getKey()[dim] >= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = node->getKey()[dim];
                result.dim = node->dim();
                invariantBroken |= result.best == node->getKey()[dim];
            }
        } else {
            // split in any other dimension.
            // First, check local key.
            double localX = node->getKey()[dim];
            if (localX >= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = localX;
                result.dim = node->dim();
                invariantBroken |= result.best == localX;
            }
            if (node->getLo() != nullptr) {
                removeMaxLeaf(node->getLo(), node, dim, result);
            }
            if (node->getHi() != nullptr) {
                removeMaxLeaf(node->getHi(), node, dim, result);
            }
        }
    }

    /**
     * Reinsert the key.
     * @param oldKey old key
     * @param newKey new key
     * @return the value associated with the key or 'null' if the key was not found.
     */
  public:
    size_t relocate(const Key& oldKey, const Key& newKey, const T& value) {
        if (root_ == nullptr) {
            return 0;
        }
        // TODO check results
        erase(oldKey, value);
        insert(newKey, value);
        return 1;
    }

    /**
     * Get the number of key-value pairs in the tree.
     * @return the size
     */
    int size() {
        return size_;
    }

    bool empty() {
        return size_ == 0;
    }

    /**
     * Removes all elements from the tree.
     */
    void clear() {
        size_ = 0;
        delete root_;
        root_ = nullptr;
        invariantBroken = false;
    }

    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(QueryBox query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        auto it = begin_query(query_box, std::forward<FILTER>(filter));
        // TODO move FILTER into begin_query
        auto end = this->end();
        while (it != end) {
            Key k = it._node()->point();
            if (filter.IsEntryValid(k, *it)) {
                callback(k, *it);
            }
            ++it;
        }
    }

    /**
     * Query the tree, returning all points in the axis-aligned rectangle between 'min' and 'max'.
     * @param min lower left corner of query
     * @param max upper right corner of query
     * @return all entries in the rectangle
     */
    template <typename FILTER = FilterNoOp>
    auto begin_query(const QueryBox& query_box, FILTER&& filter = FILTER()) const {
        return KDIterator<Key, T, FILTER>(
            root_, query_box.min(), query_box.max(), std::forward<FILTER>(filter));
    }

    auto end() const {
        return KDIterator<Key, T, FilterNoOp>();
    }

    template <typename DISTANCE, typename FILTER = FilterNoOp>
    auto begin_knn_query(
        size_t k,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        auto result = knnQuery(
            center, k, std::forward<DISTANCE>(distance_function), std::forward<FILTER>(filter));
        // TODO pass in directly w/o move()
        return KDIteratorKnn<Key, T>(std::move(result));
    }

    /**
     * 1-nearest neighbor query.
     * @param center The point for which the nearest neighbors are requested
     * @return Nearest neighbor
     */
    KDEntryDist<Key, T> nnQuery(const Key& center) {
        if (root_ == nullptr) {
            return nullptr;
        }
        KDEntryDist<Key, T> candidate = new KDEntryDist(nullptr, SCALAR_MAX);
        rangeSearch1NN(root_, center, candidate, SCALAR_MAX);
        return candidate;
    }

  private:
    double rangeSearch1NN(
        Node<Key, T>* node, const Key& center, KDEntryDist<Key, T> candidate, double maxRange) {
        dimension_t dim = node->dim();
        if (node->getLo() != nullptr &&
            (center[dim] < node->getKey()[dim] || node->getHi() == nullptr)) {
            // go down
            maxRange = rangeSearch1NN(node->getLo(), center, candidate, maxRange);
            // refine result
            if (center[dim] + maxRange >= node->getKey()[dim]) {
                maxRange = addCandidate(node, center, candidate, maxRange);
                if (node->getHi() != nullptr) {
                    maxRange = rangeSearch1NN(node->getHi(), center, candidate, maxRange);
                }
            }
        } else if (node->getHi() != nullptr) {
            // go down
            maxRange = rangeSearch1NN(node->getHi(), center, candidate, maxRange);
            // refine result
            if (center[dim] <= node->getKey()[dim] + maxRange) {
                maxRange = addCandidate(node, center, candidate, maxRange);
                if (node->getLo() != nullptr) {
                    maxRange = rangeSearch1NN(node->getLo(), center, candidate, maxRange);
                }
            }
        } else {
            // leaf -> first (probably best) match!
            maxRange = addCandidate(node, center, candidate, maxRange);
        }
        return maxRange;
    }

    double addCandidate(
        Node<Key, T>* node,
        const Key& center,
        const KDEntryDist<Key, T>& candidate,
        double maxRange) {
        double dist = distance(center, node->getKey());
        if (dist >= maxRange) {
            // don't add if too far away
            // don't add if we already have an equally good result
            return maxRange;
        }
        candidate.set(node, dist);
        return dist;
    }

  public:
    template <typename DISTANCE, typename FILTER = FilterNoOp>
    auto knnQuery(
        const Key& center,
        size_t k,
        DISTANCE&& distance_fn = DISTANCE(),
        FILTER&& filter_fn = FILTER()) const {
        if (root_ == nullptr) {
            return std::vector<KDEntryDist<Key, T>>(0);
        }
        std::vector<KDEntryDist<Key, T>> candidates{};
        candidates.reserve(k);
        DISTANCE dist_fn{std::forward<DISTANCE>(distance_fn)};
        FILTER filt_fn{std::forward<FILTER>(filter_fn)};
        rangeSearchKNN(root_, center, candidates, k, SCALAR_MAX, dist_fn, filt_fn);
        return candidates;
    }

  private:
    template <typename DISTANCE, typename FILTER = FilterNoOp>
    double rangeSearchKNN(
        Node<Key, T>* node,
        const Key& center,
        std::vector<KDEntryDist<Key, T>>& candidates,
        size_t k,
        double maxRange,
        DISTANCE& dist_fn,
        FILTER& filter) const {
        dimension_t dim = node->dim();
        if (node->getLo() != nullptr &&
            (center[dim] < node->getKey()[dim] || node->getHi() == nullptr)) {
            // go down
            maxRange =
                rangeSearchKNN(node->getLo(), center, candidates, k, maxRange, dist_fn, filter);
            // refine result
            if (center[dim] + maxRange >= node->getKey()[dim]) {
                maxRange = addCandidate(node, center, candidates, k, maxRange, dist_fn, filter);
                if (node->getHi() != nullptr) {
                    maxRange = rangeSearchKNN(
                        node->getHi(), center, candidates, k, maxRange, dist_fn, filter);
                }
            }
        } else if (node->getHi() != nullptr) {
            // go down
            maxRange =
                rangeSearchKNN(node->getHi(), center, candidates, k, maxRange, dist_fn, filter);
            // refine result
            if (center[dim] <= node->getKey()[dim] + maxRange) {
                maxRange = addCandidate(node, center, candidates, k, maxRange, dist_fn, filter);
                if (node->getLo() != nullptr) {
                    maxRange = rangeSearchKNN(
                        node->getLo(), center, candidates, k, maxRange, dist_fn, filter);
                }
            }
        } else {
            // leaf -> first (probably best) match!
            maxRange = addCandidate(node, center, candidates, k, maxRange, dist_fn, filter);
        }
        return maxRange;
    }

    template <typename DISTANCE, typename FILTER = FilterNoOp>
    double addCandidate(
        Node<Key, T>* node,
        const Key& center,
        std::vector<KDEntryDist<Key, T>>& candidates,
        size_t k,
        double maxRange,
        DISTANCE& dist_fn,
        FILTER& filter) const {
        // add ?
        double dist = dist_fn(center, node->getKey());
        if (dist > maxRange) {
            // don't add if too far away
            return maxRange;
        }
        if (dist == maxRange && candidates.size() >= k) {
            // don't add if we already have enough equally good results.
            return maxRange;
        }
        if (!filter.IsEntryValid(node->getKey(), node->getValue())) {
            return maxRange;
        }
        //        KDEntryDist<Key, T> cand{};
        if (candidates.size() >= k) {
            //            cand = candidates.erase(candidates.end() - 1);
            candidates.erase(candidates.end() - 1);
            //            cand.set(node, dist);
            //        } else {
            //            cand = new KDEntryDist(node, dist);
        }
        //        auto compKnn = [](const KDEntryDist<Key, T>& p1, const KDEntryDist<Key, T>& p2) {
        //            double deltaDist = p1.dist() - p2.dist();
        //            return deltaDist < 0 ? -1 : (deltaDist > 0 ? 1 : 0);
        //        };
        //        int insertionPos = std::binary_search(candidates.begin(), candidates.end(), cand,
        //        compKnn); insertionPos = insertionPos >= 0 ? insertionPos : -(insertionPos + 1);
        //        candidates.emplace(insertionPos, cand);
        //        return candidates.size() < k ? maxRange : candidates.get(candidates.size() -
        //        1).dist();
        // TODO upper_bound to avoid moving things?
        auto insertionPos = std::lower_bound(
            candidates.begin(),
            candidates.end(),
            dist,
            [](const KDEntryDist<Key, T>& p1, double dist) { return p1.dist() < dist; });
        // insertionPos = insertionPos >= 0 ? insertionPos : -(insertionPos + 1);
        candidates.emplace(insertionPos, node, dist);
        return candidates.size() < k ? maxRange : candidates.back().dist();
    }

  private:
    /**
     * Returns a printable list of the tree.
     * @return the tree as String
     */
  public:
    [[nodiscard]] std::string to_string() const {
        std::ostringstream os;
        toStringTree(os, root_, 0);
        return os.str();
    }

    std::ostream& operator<<(std::ostream& out) {
        if (root_ == nullptr) {
            out << "empty tree";
        } else {
            // TODO this is weird, just use << ?!?
            toStringTree(out, root_, 0);
        }
        return out;
    }

  private:
    void toStringTree(std::ostream& out, Node<Key, T>* node, int depth) const {
        std::string prefix = "";
        for (int i = 0; i < depth; i++) {
            prefix += ".";
        }
        prefix += " ";
        out << prefix << node->point();
        out << " v=" << node->value();
        out << " l/r=";
        if (node->getLo() == nullptr) {
            out << "nullptr";
        } else {
            out << node->getLo()->point();
        }
        out << "/";
        if (node->getHi() == nullptr) {
            out << "nullptr";
        } else {
            out << node->getHi()->point();
        }
        out << std::endl;
        if (node->getLo() != nullptr) {
            toStringTree(out, node->getLo(), depth + 1);
        }
        if (node->getHi() != nullptr) {
            toStringTree(out, node->getHi(), depth + 1);
        }
    }

  public:
    KDStats stats() {
        KDStats s(this);
        if (root_ != nullptr) {
            root_->stats(s, 0);
        }
        return s;
    }

    void check_consistency() const {
        if (root_ != nullptr) {
            Key lower_bound;
            Key upper_bound;
            for (size_t d = 0; d < lower_bound.size(); ++d) {
                lower_bound[d] = SCALAR_MIN;
                upper_bound[d] = SCALAR_MAX;
            }
            size_t n = 0;
            root_->check_consistency(lower_bound, upper_bound, n, lower_bound.size() - 1);
            assert(n == size_);
        }
    }

    size_t dims() {
        return dims_;
    }

    auto begin() {
        return begin_query(
            {{SCALAR_MIN, SCALAR_MIN, SCALAR_MIN}, {SCALAR_MAX, SCALAR_MAX, SCALAR_MAX}});
    }

    auto query1NN(const Key& center) {
        return nnQuery(center);
    }

  private:
    const size_t dims_;
    size_t size_ = 0;
    // During insertion, the tree maintains an invariant that if two points have the
    // same value in any dimension, then one key is never in the 'lower' branch of the other.
    // This allows very efficient look-up because we have to follow only a single path.
    // Unfortunately, removing keys (and moving up keys in the tree) may break this invariant.
    //
    // The cost of breaking the invariant is that there may be two branches that contain the same
    // key. This makes searches more expensive simply because the code gets more complex.
    //
    // One solution to maintain the invariant is to repair the invariant by moving all points
    // with the same value into the 'upper' part of the point that was moved up.
    // Identifying these points is relatively cheap, we can do this while searching for
    // the min/max during removal. However, _moving_ these point may be very expensive.
    //
    // Our (pragmatic) solution is to identify situations when the invariant gets broken.
    // When it is not broken, we use the simple search. If it gets broken, we use the slower search.
    // This is especially useful in scenarios where 'remove()' is not required or where
    // points have never the same values (such as for physical measurements or other experimental
    // results).
    bool invariantBroken = false;

    Node<Key, T>* root_ = nullptr;
};

}  // namespace tinspin

#endif