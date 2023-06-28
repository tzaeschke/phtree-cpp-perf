// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef TINSPIN_KD_TREE_H
#define TINSPIN_KD_TREE_H

#include "include/phtree/common/bpt_priority_queue.h"
#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include "min_max_vector_heap.h"
#include "src/util/ph-util.h"
#include <iostream>
#include <queue>
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
bool isEnclosed(const Key& point, const Key& min, const Key& max) {
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
    KDEntryDist() : entry_{nullptr}, distance_{0} {};
    KDEntryDist(Node<Key, T>* node, double dist) : entry_{node}, distance_{dist} {}

    void set(Node<Key, T>* node, double dist) {
        entry_ = node;
        distance_ = dist;
    }

    [[nodiscard]] double dist() const noexcept {
        return distance_;
    }

    const Key& key() const noexcept {
        return entry_->key();
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
    Node(const Key& key, const T& value, size_t split_dim)
    : coordinate_{key}, value_{value}, split_dim_{split_dim} {}

    ~Node() {
        delete left_;
        delete right_;
    }

    Node* get_closest_node_or_add_entry(const Key& p, const T& value) {
        // Find best sub-node.
        // If there is no node, we create one and return null
        if (p[split_dim_] >= coordinate_[split_dim_]) {
            if (right_ != nullptr) {
                return right_;
            }
            right_ = new Node(p, value, (split_dim_ + 1) % p.size());
            return nullptr;
        }
        if (left_ != nullptr) {
            return left_;
        }
        left_ = new Node(p, value, (split_dim_ + 1) % p.size());
        return nullptr;
    }

    const Key& key() const noexcept {  // backport: -> remove duplicate point/key methods
        return coordinate_;
    }

    const T& value() const noexcept {
        return value_;
    }

    Node* left() const noexcept {
        return left_;
    }

    Node* right() const noexcept {
        return right_;
    }

    void set_left(Node* left) noexcept {
        delete left_;  // deleting null is fine.
        left_ = left;
    }

    void set_right(Node* right) noexcept {
        delete right_;  // deleting null is fine.
        right_ = right;
    }

    void set(const Key& key, const T& value) noexcept {
        coordinate_ = key;
        value_ = value;
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
        assert(split_dim_ == (d_parent + 1u) % coordinate_.size());
        for (size_t d = 0; d < lower_bound.size(); ++d) {
            assert(coordinate_[d] >= lower_bound[d]);
            assert(coordinate_[d] <= upper_bound[d]);
        }

        ++n;
        if (left_ != nullptr) {
            Key upper_bound_2{upper_bound};
            upper_bound_2[split_dim_] =
                std::min(upper_bound_2[split_dim_], coordinate_[split_dim_]);
            left_->check_consistency(lower_bound, upper_bound_2, n, split_dim_);
        }
        if (right_ != nullptr) {
            Key lower_bound_2{lower_bound};
            lower_bound_2[split_dim_] =
                std::max(lower_bound_2[split_dim_], coordinate_[split_dim_]);
            right_->check_consistency(lower_bound_2, upper_bound, n, split_dim_);
        }
    }

    [[nodiscard]] bool is_leaf() const noexcept {
        return left_ == nullptr && right_ == nullptr;
    }

    [[nodiscard]] size_t dim() const noexcept {
        return split_dim_;
    }

  private:
    Key coordinate_;
    T value_;
    Node* left_ = nullptr;
    Node* right_ = nullptr;
    size_t split_dim_;
};

template <typename Key, typename T>
class KDIteratorBase {
  public:
    explicit KDIteratorBase() noexcept : node_{nullptr} {}

    inline auto& operator*() const noexcept {
        assert(node_ != nullptr);
        return node_->value();
    }

    inline auto* operator->() const noexcept {
        assert(node_ != nullptr);
        return &node_->value();
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
    [[nodiscard]] bool IsEnd() const noexcept {
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
    using depth_t = std::uint32_t;
    struct IteratorPos {
        Node<Key, T>* node_;
        depth_t depth_;
        bool doLeft, doKey, doRight;

        void set(Node<Key, T>* node, const Key& min, const Key& max, depth_t depth) {
            node_ = node;
            depth_ = depth;
            const Key& key = node_->key();
            auto dims = min.size();
            dimension_t dim = depth % dims;
            doLeft = min[dim] <= key[dim];
            doRight = max[dim] >= key[dim];
            doKey = true;
        }
    };

  public:
    // end()
    KDIterator() : KDIteratorBase<Key, T>(), stack_{}, min_{}, max_{} {}

    // begin_query()
    template <typename F = FilterNoOp>
    KDIterator(Node<Key, T>* root, const Key& min, const Key& max, F&& filter = F())
    : KDIteratorBase<Key, T>(), stack_{}, min_{min}, max_{max}, filter_(std::forward<F>(filter)) {
        if (root != nullptr) {
            prepareAndPush(root, min, max, 0);
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
        while (!stack_.empty()) {
            IteratorPos& itPos = stack_.top();
            Node<Key, T>* node = itPos.node_;
            if (itPos.doLeft && node->left() != nullptr) {
                itPos.doLeft = false;
                prepareAndPush(node->left(), min_, max_, itPos.depth_ + 1);
                continue;
            }
            if (itPos.doKey) {
                itPos.doKey = false;
                if (isEnclosed(node->key(), min_, max_) &&
                    filter_.IsEntryValid(node->key(), node->value())) {
                    this->SetCurrentResult(node);
                    return;
                }
            }
            if (itPos.doRight && node->right() != nullptr) {
                itPos.doRight = false;
                prepareAndPush(node->right(), min_, max_, itPos.depth_ + 1);
                continue;
            }
            stack_.pop();
        }
        this->SetFinished();
    }

    IteratorPos prepareAndPush(Node<Key, T>* node, const Key& min, const Key& max, depth_t depth) {
        auto& n = stack_.push();
        n.set(node, min, max, depth);
        return n;
    }

  private:
    tinspin::SimpleStack<IteratorPos> stack_;
    Key min_;
    Key max_;
    FILTER filter_;
};

template <typename Key, typename T>
class KDIteratorKnn : public KDIteratorBase<Key, T> {
    using Candidates = std::vector<KDEntryDist<Key, T>>;
    using CandidatesIter = decltype(Candidates{}.begin());

  public:
    explicit KDIteratorKnn(Candidates&& result) noexcept
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
        return iter_->key();
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

template <typename Key, typename Value>
struct NodeDist {
    NodeDist(const Node<Key, Value>* node, const Key& min, const Key& max)
    : dist_{0}, node_{node}, min_{min}, max_{max} {}

    NodeDist(double dist, const Node<Key, Value>* node, const Key& min, const Key& max)
    : dist_{dist}, node_{node}, min_{min}, max_{max} {}

    double dist_;
    const Node<Key, Value>* node_{nullptr};
    Key min_;
    Key max_;
};

template <typename Key, typename Value>
struct ValueDist {
    ValueDist(double distance, const Node<Key, Value>* node) noexcept
    : dist_{distance}, node_{node} {}
    double dist_;
    const Node<Key, Value>* node_{nullptr};
};

template <typename Key, typename Value, typename DISTANCE, typename FILTER>
class KDIteratorKnnHS : public KDIteratorBase<Key, Value> {
    using SCALAR = std::remove_reference_t<decltype(Key{}[0])>;
    using EntryT = Node<Key, Value>;
    using NodeDistT = NodeDist<Key, Value>;
    using ValueDistT = ValueDist<Key, Value>;

    struct CompareNodeDistByDistance {
        bool operator()(const NodeDistT& left, const NodeDistT& right) const {
            return left.dist_ > right.dist_;
        };
    };
    struct CompareValueDistByDistance {
        bool operator()(const ValueDistT& left, const ValueDistT& right) const {
            return left.dist_ > right.dist_;
        };
    };

  public:
    template <typename DIST, typename F>
    explicit KDIteratorKnnHS(
        const EntryT* root, const Key& center, size_t min_results, DIST&& dist, F&& filter)
    : KDIteratorBase<Key, Value>()
    , center_{center}
    , current_distance_{std::numeric_limits<double>::max()}
    , remaining_{min_results}
    , filter_{std::forward<F>(filter)}
    , distance_(std::forward<DIST>(dist)) {
        if (min_results <= 0 || root == nullptr) {
            this->SetFinished();
            return;
        }

        // Initialize queue, use d=0 because every imaginable point lies inside the root Node
        Key min{};
        Key max{};
        for (dimension_t d = 0; d < min.size(); ++d) {
            min[d] = -std::numeric_limits<SCALAR>::infinity();
            max[d] = std::numeric_limits<SCALAR>::infinity();
        }
        // queue_n_.emplace(NodeDistT{0, root});
        queue_n_.emplace(root, min, max);
        FindNextElement();
    }

    [[nodiscard]] double distance() const {
        return current_distance_;
    }

    const Key& first() {
        return this->_node()->key();
    }

    KDIteratorKnnHS& operator++() noexcept {
        FindNextElement();
        return *this;
    }

    KDIteratorKnnHS operator++(int) noexcept {
        KDIteratorKnnHS iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void FindNextElement() {
        while (remaining_ > 0 && !(queue_n_.empty() && queue_v_.empty())) {
            bool use_v = !queue_v_.empty();
            if (use_v && !queue_n_.empty()) {
                use_v = queue_v_.top().dist_ <= queue_n_.top().dist_;
            }
            if (use_v) {
                // data entry
                auto& result = queue_v_.top();
                --remaining_;
                this->SetCurrentResult(const_cast<Node<Key, Value>*>(result.node_));  // TODO cast?
                current_distance_ = result.dist_;
                queue_v_.pop();
                return;
            } else {
                // inner node
                auto& top = queue_n_.top();  // This creates a copy!
                auto& node = *top.node_;
                auto d_node = top.dist_;
                auto min = top.min_;  // copy
                auto max = top.max_;  // copy
                queue_n_.pop();

                if (d_node > max_node_dist_ && queue_v_.size() >= remaining_) {
                    // ignore this node
                    continue;
                }

                if (filter_.IsEntryValid(node.key(), node.value())) {
                    double d = distance_(center_, node.key());
                    // Using '<=' allows dealing with infinite distances.
                    if (d <= max_node_dist_) {
                        queue_v_.emplace(d, &node);
                        if (queue_v_.size() >= remaining_) {
                            if (queue_v_.size() > remaining_) {
                                queue_v_.pop_max();
                            }
                            double d_max = queue_v_.top_max().dist_;
                            max_node_dist_ = std::min(max_node_dist_, d_max);
                        }
                    }
                }
                // left
                if (node.left() != nullptr) {
                    // create_entry(top, false, min, max);
                    create_entry_left(&node, min, max);
                }
                // right
                if (node.right() != nullptr) {
                    // create_entry(top, true, min, max);
                    create_entry_right(&node, min, max);
                }
            }
        }
        this->SetFinished();
        current_distance_ = std::numeric_limits<double>::max();
        // TODO test all with 20D/30D
    }

    void create_entry_right(const Node<Key, Value>* node, const Key& min, const Key& max) noexcept {
        auto split_dim = node->dim();
        assert(node->key()[split_dim] <= max[split_dim]);
        assert(node->key()[split_dim] >= min[split_dim]);
        auto min_new = min;  // copy
        min_new[split_dim] = node->key()[split_dim];
        auto dist = DistanceToNode(min_new, max);
        if (dist <= max_node_dist_) {
            queue_n_.emplace(dist, node->right(), std::move(min_new), max);
        }
    }

    void create_entry_left(const Node<Key, Value>* node, const Key& min, const Key& max) noexcept {
        auto split_dim = node->dim();
        assert(node->key()[split_dim] <= max[split_dim]);
        assert(node->key()[split_dim] >= min[split_dim]);
        auto max_new = max;  // copy
        max_new[split_dim] = node->key()[split_dim];
        auto dist = DistanceToNode(min, max_new);
        if (dist <= max_node_dist_) {
            queue_n_.emplace(dist, node->left(), min, std::move(max_new));
        }
    }

    double DistanceToNode(const Key& node_min, const Key& node_max) {
        Key buf;
        // The following calculates the point inside the node that is closest to center_.
        for (dimension_t i = 0; i < node_min.size(); ++i) {
            // if center_[i] is outside the node, return distance to the closest edge,
            // otherwise return center_[i] itself (assume possible distance=0)
            SCALAR min = node_min[i];
            SCALAR max = node_max[i];
            buf[i] = min > center_[i] ? min : (max < center_[i] ? max : center_[i]);
        }
        return distance_(center_, buf);
    }

  private:
    const Key center_;
    double current_distance_;
    size_t remaining_;
    // std::priority_queue<NodeDistT, std::vector<NodeDistT>, CompareNodeDistByDistance> queue_n_;
    aux::min_max_vector_heap<NodeDistT, CompareNodeDistByDistance> queue_n_;
    aux::min_max_vector_heap<ValueDistT, CompareValueDistByDistance> queue_v_;
    //    ::phtree::bptree::detail::priority_queue<EntryDistT, CompareEntryDistByDistance> queue_n_;
    //    ::phtree::bptree::detail::priority_queue<EntryDistT, CompareEntryDistByDistance> queue_v_;
    FILTER filter_;
    DISTANCE distance_;
    double max_node_dist_ = std::numeric_limits<double>::infinity();
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
        Node<Key, T>* node{nullptr};
        Node<Key, T>* nodeParent{nullptr};
        double best{0};
        dimension_t dim{0};
    };

  public:
    using KeyInternal = Key;

    KDTree() {
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
        while ((n = n->get_closest_node_or_add_entry(key, value)) != nullptr)
            ;
    }

    void insert(const Key& key, const T& value) {
        size_++;
        if (root_ == nullptr) {
            root_ = new Node(key, value, 0);
            return;
        }
        Node<Key, T>* n = root_;
        while ((n = n->get_closest_node_or_add_entry(key, value)) != nullptr)
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
            const Key& nodeKey = n->key();
            double nodeX = nodeKey[n->dim()];
            double keyX = key[n->dim()];
            if (keyX == nodeX && key == nodeKey && pred_fn(n->value())) {
                resultDepth.dim = n->dim();
                resultDepth.nodeParent = parent;
                return n;
            }
            parent = n;
            n = (keyX >= nodeX) ? n->right() : n->left();
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
            auto& nodeKey = n->key();
            double nodeX = nodeKey[n->dim()];
            double keyX = key[n->dim()];
            if (keyX == nodeX) {
                if (key == nodeKey && pred_fn(n->value())) {
                    resultDepth.dim = n->dim();
                    resultDepth.nodeParent = parent;
                    return n;
                }
                // Broken invariant? We need to check the 'lower' part as well...
                if (n->left() != nullptr) {
                    Node<Key, T>* n2 = findNodeExactSlow(key, n->left(), n, resultDepth, pred_fn);
                    if (n2 != nullptr) {
                        return n2;
                    }
                }
            }
            parent = n;
            n = (keyX >= nodeX) ? n->right() : n->left();
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

        invariantBroken = true;

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
            return 1;
        }

        // find replacement
        while (eToRemove != nullptr && !eToRemove->is_leaf()) {
            // recurse
            auto dim = removeResult.dim;
            removeResult.node = nullptr;
            if (eToRemove->right() != nullptr) {
                // get replacement from right
                // This is preferable, because it cannot break the invariant
                removeResult.best = SCALAR_MAX;
                removeMinLeaf(eToRemove->right(), eToRemove, dim, removeResult);
            } else if (eToRemove->left() != nullptr) {
                // get replacement from left
                removeResult.best = SCALAR_MIN;
                removeMaxLeaf(eToRemove->left(), eToRemove, dim, removeResult);
            }
            eToRemove->set(removeResult.node->key(), removeResult.node->value());
            eToRemove = removeResult.node;
        }
        // leaf node
        Node<Key, T>* parent = removeResult.nodeParent;
        if (parent != nullptr) {
            if (parent->left() == eToRemove) {
                parent->set_left(nullptr);
            } else if (parent->right() == eToRemove) {
                parent->set_right(nullptr);
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
            if (node->left() != nullptr) {
                removeMinLeaf(node->left(), node, dim, result);
            } else if (node->key()[dim] <= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = node->key()[dim];
                result.dim = node->dim();
            }
        } else {
            // split in any other dimension.
            // First, check local key.
            double localX = node->key()[dim];
            if (localX <= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = localX;
                result.dim = node->dim();
            }
            if (node->left() != nullptr) {
                removeMinLeaf(node->left(), node, dim, result);
            }
            if (node->right() != nullptr) {
                removeMinLeaf(node->right(), node, dim, result);
            }
        }
    }

    void removeMaxLeaf(
        Node<Key, T>* node, Node<Key, T>* parent, dimension_t dim, RemoveResult& result) {
        // Split in 'interesting' dimension
        if (dim == node->dim()) {
            // We strictly look for leaf nodes with left==null
            if (node->right() != nullptr) {
                removeMaxLeaf(node->right(), node, dim, result);
            } else if (node->key()[dim] >= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = node->key()[dim];
                result.dim = node->dim();
            }
        } else {
            // split in any other dimension.
            // First, check local key.
            double localX = node->key()[dim];
            if (localX >= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = localX;
                result.dim = node->dim();
            }
            if (node->left() != nullptr) {
                removeMaxLeaf(node->left(), node, dim, result);
            }
            if (node->right() != nullptr) {
                removeMaxLeaf(node->right(), node, dim, result);
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
    size_t size() {
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
            Key k = it._node()->key();
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
    auto begin_knn_query1(
        size_t k,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        auto result = knnQuery(
            center, k, std::forward<DISTANCE>(distance_function), std::forward<FILTER>(filter));
        // TODO pass in directly w/o move()
        return KDIteratorKnn<Key, T>(std::move(result));
    }

    template <typename DISTANCE, typename FILTER = FilterNoOp>
    auto begin_knn_query(
        size_t k,
        const Key& center,
        DISTANCE&& distance_function = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        return KDIteratorKnnHS<Key, T, DISTANCE, FILTER>(
            root_,
            center,
            k,
            std::forward<DISTANCE>(distance_function),
            std::forward<FILTER>(filter));
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
        if (node->left() != nullptr &&
            (center[dim] < node->key()[dim] || node->right() == nullptr)) {
            // go down
            maxRange = rangeSearch1NN(node->left(), center, candidate, maxRange);
            // refine result
            if (center[dim] + maxRange >= node->key()[dim]) {
                maxRange = addCandidate(node, center, candidate, maxRange);
                if (node->right() != nullptr) {
                    maxRange = rangeSearch1NN(node->right(), center, candidate, maxRange);
                }
            }
        } else if (node->right() != nullptr) {
            // go down
            maxRange = rangeSearch1NN(node->right(), center, candidate, maxRange);
            // refine result
            if (center[dim] <= node->key()[dim] + maxRange) {
                maxRange = addCandidate(node, center, candidate, maxRange);
                if (node->left() != nullptr) {
                    maxRange = rangeSearch1NN(node->left(), center, candidate, maxRange);
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
        double dist = distance(center, node->key());
        if (dist >= maxRange) {
            // don't add if too far away
            // don't add if we already have an equally good result
            return maxRange;
        }
        candidate.set(node, dist);
        return dist;
    }

    using CandidateContainerT = std::vector<KDEntryDist<Key, T>>;  // TODO use below...?
    // using CandidateContainerT = std::map<KDEntryDist<Key, T>>;
    // using CandidateContainerT = std::priority_queue<KDEntryDist<Key, T>>;

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
        if (node->left() != nullptr &&
            (center[dim] < node->key()[dim] || node->right() == nullptr)) {
            // go down
            maxRange =
                rangeSearchKNN(node->left(), center, candidates, k, maxRange, dist_fn, filter);
            // refine result
            if (center[dim] + maxRange >= node->key()[dim]) {
                maxRange = addCandidate(node, center, candidates, k, maxRange, dist_fn, filter);
                if (node->right() != nullptr) {
                    maxRange = rangeSearchKNN(
                        node->right(), center, candidates, k, maxRange, dist_fn, filter);
                }
            }
        } else if (node->right() != nullptr) {
            // go down
            maxRange =
                rangeSearchKNN(node->right(), center, candidates, k, maxRange, dist_fn, filter);
            // refine result
            if (center[dim] <= node->key()[dim] + maxRange) {
                maxRange = addCandidate(node, center, candidates, k, maxRange, dist_fn, filter);
                if (node->left() != nullptr) {
                    maxRange = rangeSearchKNN(
                        node->left(), center, candidates, k, maxRange, dist_fn, filter);
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
        double dist = dist_fn(center, node->key());
        if (dist > maxRange) {
            // don't add if too far away
            return maxRange;
        }
        if (dist == maxRange && candidates.size() >= k) {
            // don't add if we already have enough equally good results.
            return maxRange;
        }
        if (!filter.IsEntryValid(node->key(), node->value())) {
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
        out << prefix << node->key();
        out << " v=" << node->value();
        out << " l/r=";
        if (node->left() == nullptr) {
            out << "nullptr";
        } else {
            out << node->left()->key();
        }
        out << "/";
        if (node->right() == nullptr) {
            out << "nullptr";
        } else {
            out << node->right()->key();
        }
        out << std::endl;
        if (node->left() != nullptr) {
            toStringTree(out, node->left(), depth + 1);
        }
        if (node->right() != nullptr) {
            toStringTree(out, node->right(), depth + 1);
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

    auto begin() {
        return begin_query(
            {{SCALAR_MIN, SCALAR_MIN, SCALAR_MIN}, {SCALAR_MAX, SCALAR_MAX, SCALAR_MAX}});
    }

    auto query1NN(const Key& center) {
        return nnQuery(center);
    }

  private:
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