// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef TINSPIN_KD_TREE_H
#define TINSPIN_KD_TREE_H

#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include <iostream>
#include <unordered_set>

namespace tinspin {

using namespace improbable::phtree;

namespace {
template<typename Key>
static bool isEnclosed(const Key& point, const Key& min, const Key& max) {
    for (int i = 0; i < point.length; i++) {
        if (point[i] < min[i] || point[i] > max[i]) {
            return false;
        }
    }
    return true;
}
}


template <typename Key, typename T>
class Node;

template <typename T, typename SCALAR = double>
class KDTree;

struct KDStats {
    size_t nNodes = 0;
    size_t maxDepth = 0;
};

template <typename Key, typename T>
class KDEntryDist {  // implements PointEntryDist<T> {

  private:
    Node<Key, T>* entry_;

  private:
    double distance_;

  public:
    KDEntryDist(Node<Key, T>* node, double dist) {
        entry_ = node;
        distance_ = dist;
    }

    void set(Node<Key, T>* node, double dist) {
        entry_ = node;
        distance_ = dist;
    }

  public:
    double dist() {
        return distance_;
    }

  public:
    static const QEntryComparator COMP = new QEntryComparator();

    static class QEntryComparator implements Comparator < KDEntryDist < ? >> {
        /**
         * Compares the two specified MBRs according to
         * the sorting dimension and the sorting co-ordinate for the dimension
         * of this Comparator.
         *
         * @param o1 the first SpatialPoint
         * @param o2 the second SpatialPoint
         * @return a negative integer, zero, or a positive integer as the
         *         first argument is less than, equal to, or greater than the
         *         second.
         */
      public: int compare(KDEntryDist<?> o1, KDEntryDist<?> o2) {
            double d = o1.dist() - o2.dist();
            return d < 0 ? -1 : (d > 0 ? 1 : 0);
        }
    }

    public : const Key&
             point() {
        return entry_->point();
    }

  public:
    const T& value() {
        return entry_->value();
    }
};

/**
 * Node class for the kdtree.
 */
template <typename Coord, typename T>
class Node {  //: PointEntry<T> {
  private:
    Coord coordinate_;

  private:
    T value_;
    Node* left_ = nullptr;
    Node* right_ = nullptr;
    int dim_;

  public:
    Node(const Coord& p, const T& value, int dim) {
        coordinate_ = p;
        value_ = value;
        dim_ = dim;
    }

    Node* getClosestNodeOrAddPoint(const Coord& p, const T& value, int dims) const {
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

    const Coord& getKey() const noexcept {
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
        left_ = left;
    }

    void setRight(Node* right) noexcept {
        right_ = right;
    }

    void setKeyValue(const Coord& key, const T& value) noexcept {
        coordinate_ = key;
        value_ = value;
    }

    const Coord& point() const noexcept {
        return coordinate_;
    }

    const T& value() const noexcept {
        return value_;
    }

    void checkNode(KDStats& s, int depth) {
        s.nNodes++;
        if (depth > s.maxDepth) {
            s.maxDepth = depth;
        }
        if (left_ != nullptr) {
            left_->checkNode(s, depth + 1);
        }
        if (right_ != nullptr) {
            right_->checkNode(s, depth + 1);
        }
    }

    //    char* toString() {
    //        return "center=" + Arrays.toString(point()) + " " + (void*)(this);
    //    }

    bool isLeaf() const noexcept {
        return left_ == nullptr && right_ == nullptr;
    }

    int getDim() const noexcept {
        return dim_;
    }
};

/**
 * Resetable query iterator.
 *
 * @param <T> Value type
 */
template <typename Key, typename T>
class KDIterator {  //}; implements QueryIterator<PointEntry<T>> {

    class IteratorPos {
        Node<Key, T>* node_;
        int depth_;
        bool doLeft, doKey, doRight;

        void set(Node<Key, T>* node, const Key& min, const Key& max, int depth, int dims) {
            node_ = node;
            depth_ = depth;
            const Key& key = node_->getKey();
            int pos = depth % dims;
            doLeft = min[pos] < key[pos];
            doRight = max[pos] > key[pos];
            doKey = doLeft || doRight || key[pos] == min[pos] || key[pos] == max[pos];
        }
    };

  private:
    class IteratorStack {
      private:
        const std::vector<IteratorPos> stack{};

      private:
        int size = 0;

        IteratorStack() {}

        bool isEmpty() {
            return size == 0;
        }

        IteratorPos prepareAndPush(
            Node<Key, T>* node, const Key& min, const Key& max, int depth, int dims) {
            if (size == stack.size()) {
                stack.emplace_back();
            }
            IteratorPos& ni = stack.get(size++);

            ni.set(node, min, max, depth, dims);
            return ni;
        }

        IteratorPos& peek() {
            return stack.get(size - 1);
        }

        IteratorPos& pop() {
            return stack.get(--size);
        }

      public:
        void clear() {
            size = 0;
        }
    };

    const KDTree<T>* tree_;
    IteratorStack stack_;
    Node<Key, T>* next_ = nullptr;
    const Key& min_;
    const Key& max_;

    KDIterator(KDTree<T>* tree, const Key& min, const Key& max) : tree_{tree}, stack_{} {
        reset(min, max);
    }

  private:
    void findNext() {
        while (!stack_.isEmpty()) {
            IteratorPos itPos = stack_.peek();
            Node<Key, T>* node = itPos.node;
            if (itPos.doLeft && node->getLo() != nullptr) {
                itPos.doLeft = false;
                stack_.prepareAndPush(node->getLo(), min_, max_, itPos.depth_ + 1, tree_->getDims());
                continue;
            }
            if (itPos.doKey) {
                itPos.doKey = false;
                if (KDTree.isEnclosed(node->getKey(), min, max)) {
                    next_ = node;
                    return;
                }
            }
            if (itPos.doRight && node->getHi() != nullptr) {
                itPos.doRight = false;
                stack_.prepareAndPush(node->getHi(), min_, max_, itPos.depth_ + 1, tree_->getDims());
                continue;
            }
            stack_.pop();
        }
        next_ = nullptr;
    }

  public:
    bool hasNext() {
        return next_ != nullptr;
    }

  public:
    Node<Key, T>* next() {
        assert(hasNext());
        Node<Key, T>* ret = next_;
        findNext();
        return ret;
    }

    /**
     * Reset the iterator. This iterator can be reused in order to reduce load on the
     * garbage collector.
     * @param min lower left corner of query
     * @param max upper right corner of query
     */

  public:
    void reset(const Key& min, const Key& max) {
        stack_.clear();
        min_ = min;
        max_ = max;
        next_ = nullptr;
        if (tree_->getRoot() != nullptr) {
            stack_.prepareAndPush(tree_->getRoot(), min, max, 0, tree_->getDims());
            findNext();
        }
    }
};

/**
 * A simple KD-Tree implementation.
 *
 * @author T. Zäschke
 *
 * @param <T> Value type
 */
template <typename T, typename SCALAR>
class KDTree {  //: PointIndex<T> {

    // TODO we could use integers but would need to fix use of infinity.
    static_assert(std::is_floating_point_v<SCALAR>);

    using Key = PhPoint<3, SCALAR>;

  public:
    static const bool DEBUG = false;

  private:
    const int dims_;
    int size_ = 0;
    int modCount_ = 0;
    long nDist1NN = 0;
    long nDistKNN = 0;
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
    // Our (pragmatic) solution is identify when the invariant gets broken.
    // When it is not broken, we use the simple search. If it gets broken, we use the slower search.
    // This is especially useful in scenarios where 'remove()' is not required or where
    // points have never the same values (such as for physical measurements or other experimental
    // results).
    bool invariantBroken = false;

    Node<Key, T>* root = nullptr;

    const PointDistanceFunction dist;

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

  private:
    KDTree(int dims, PointDistanceFunction dist) {
        if (DEBUG) {
            std::cout << "Warning: DEBUG enabled" << std::endl;
        }
        this.dims_ = dims;
        this.dist = dist != nullptr ? dist : PointDistanceFunction.L2;
    }

  public:
    static<T> KDTree<T> create(int dims) {
        return new KDTree<>(dims, PointDistanceFunction.L2);
    }

  public:
    static<T> KDTree<T> create(int dims, PointDistanceFunction dist) {
        return new KDTree<>(dims, dist);
    }

    /**
     * Insert a key-value pair.
     * @param key the key
     * @param value the value
     */
    void insert(const Key& key, const T& value) {
        size_++;
        modCount_++;
        if (root == null) {
            root = new Node<>(key, value, 0);
            return;
        }
        Node<Key, T>* n = root;
        while ((n = n.getClosestNodeOrAddPoint(key, value, dims_)) != null)
            ;
    }

    /**
     * Check whether a given key exists.
     * @param key the key to check
     * @return true iff the key exists
     */
    bool containsExact(const Key& key) {
        return findNodeExcat(key, new RemoveResult<>()) != null;
    }

    /**
     * Get the value associates with the key.
     * @param key the key to look up
     * @return the value for the key or 'null' if the key was not found
     */
  public
    const T& queryExact(const Key& key) {
        Node<Key, T>* e = findNodeExcat(key, new RemoveResult<>());
        return e == null ? null : e.getValue();
    }

  private:
    Node<Key, T>* findNodeExcat(const Key& key, RemoveResult<T> resultDepth) {
        if (root == null) {
            return null;
        }
        return invariantBroken ? findNodeExactSlow(key, root, null, resultDepth)
                               : findNodeExcatFast(key, null, resultDepth);
    }

  private:
    Node<Key, T>* findNodeExcatFast(
        const Key& key, Node<Key, T>* parent, RemoveResult<T> resultDepth) {
        Node<Key, T>* n = root;
        do {
            const Key& nodeKey = n.getKey();
            double nodeX = nodeKey[n.getDim()];
            double keyX = key[n.getDim()];
            if (keyX == nodeX && Arrays.equals(key, nodeKey)) {
                resultDepth.pos = n.getDim();
                resultDepth.nodeParent = parent;
                return n;
            }
            parent = n;
            n = (keyX >= nodeX) ? n.getHi() : n.getLo();
        } while (n != null);
        return n;
    }

  private:
    Node<Key, T>* findNodeExactSlow(
        const Key& key, Node<Key, T>* n, Node<Key, T>* parent, RemoveResult<T> resultDepth) {
        do {
            auto& nodeKey = n.getKey();
            double nodeX = nodeKey[n.getDim()];
            double keyX = key[n.getDim()];
            if (keyX == nodeX) {
                if (Arrays.equals(key, nodeKey)) {
                    resultDepth.pos = n.getDim();
                    resultDepth.nodeParent = parent;
                    return n;
                }
                // Broken invariant? We need to check the 'lower' part as well...
                if (n.getLo() != null) {
                    Node<Key, T>* n2 = findNodeExactSlow(key, n.getLo(), parent, resultDepth);
                    if (n2 != null) {
                        return n2;
                    }
                }
            }
            parent = n;
            n = (keyX >= nodeX) ? n.getHi() : n.getLo();
        } while (n != null);
        return n;
    }

    /**
     * Remove a key.
     * @param key key to remove
     * @return the value associated with the key or 'null' if the key was not found
     */
  public:
    const T& remove(const Key& key) {
        if (root == null) {
            return null;
        }

        // find
        RemoveResult<T> removeResult = new RemoveResult<>();
        Node<Key, T>* eToRemove = findNodeExcat(key, removeResult);
        if (eToRemove == null) {
            return null;
        }

        // remove
        modCount_++;
        const T& value = eToRemove.getValue();
        if (eToRemove == root && size_ == 1) {
            root = null;
            size_ = 0;
            invariantBroken = false;
        }

        // find replacement
        removeResult.nodeParent = null;
        while (eToRemove != null && !eToRemove.isLeaf()) {
            // recurse
            int pos = removeResult.pos;
            removeResult.node = null;
            // randomize search direction (modCount_)
            //			if (((modCount_ & 0x1) == 0 || eToRemove.getHi() == null) &&
            // eToRemove.getLo()
            //!= null) {
            //				//get replacement from left
            //				removeResult.best = Double.NEGATIVE_INFINITY;
            //				removeMaxLeaf(eToRemove.getLo(), eToRemove, pos, removeResult);
            //			} else if (eToRemove.getHi() != null) {
            //				//get replacement from right
            //				removeResult.best = Double.POSITIVE_INFINITY;
            //				removeMinLeaf(eToRemove.getHi(), eToRemove, pos, removeResult);
            //			}
            if (eToRemove.getHi() != null) {
                // get replacement from right
                // This is preferable, because it cannot break the invariant
                removeResult.best = Double.POSITIVE_INFINITY;
                removeMinLeaf(eToRemove.getHi(), eToRemove, pos, removeResult);
            } else if (eToRemove.getLo() != null) {
                // get replacement from left
                removeResult.best = Double.NEGATIVE_INFINITY;
                removeMaxLeaf(eToRemove.getLo(), eToRemove, pos, removeResult);
            }
            eToRemove.setKeyValue(removeResult.node->getKey(), removeResult.node->getValue());
            eToRemove = removeResult.node;
        }
        // leaf node
        Node<Key, T>* parent = removeResult.nodeParent;
        if (parent != null) {
            if (parent.getLo() == eToRemove) {
                parent.setLeft(null);
            } else if (parent.getHi() == eToRemove) {
                parent.setRight(null);
            } else {
                throw new IllegalStateException();
            }
        }
        size_--;
        return value;
    }

  private:
    static class RemoveResult<T> {
        Node<Key, T>* node = null;
        Node<Key, T>* nodeParent = null;
        double best;
        int pos;
    }

    private
    : void
      removeMinLeaf(Node<Key, T>* node, Node<Key, T>* parent, int pos, RemoveResult<T> result) {
        // Split in 'interesting' dimension
        if (pos == node->getDim()) {
            // We strictly look for leaf nodes with left==null
            //  -> left!=null means the left child is at least as small as the current node
            if (node->getLo() != null) {
                removeMinLeaf(node->getLo(), node, pos, result);
            } else if (node->getKey()[pos] <= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = node->getKey()[pos];
                result.pos = node->getDim();
            }
        } else {
            // split in any other dimension.
            // First, check local key.
            double localX = node->getKey()[pos];
            if (localX <= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = localX;
                result.pos = node->getDim();
            }
            if (node->getLo() != null) {
                removeMinLeaf(node->getLo(), node, pos, result);
            }
            if (node->getHi() != null) {
                removeMinLeaf(node->getHi(), node, pos, result);
            }
        }
    }

  private:
    void removeMaxLeaf(Node<Key, T>* node, Node<Key, T>* parent, int pos, RemoveResult<T> result) {
        // Split in 'interesting' dimension
        if (pos == node->getDim()) {
            // We strictly look for leaf nodes with left==null
            if (node->getHi() != null) {
                removeMaxLeaf(node->getHi(), node, pos, result);
            } else if (node->getKey()[pos] >= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = node->getKey()[pos];
                result.pos = node->getDim();
                invariantBroken |= result.best == node->getKey()[pos];
            }
        } else {
            // split in any other dimension.
            // First, check local key.
            double localX = node->getKey()[pos];
            if (localX >= result.best) {
                result.node = node;
                result.nodeParent = parent;
                result.best = localX;
                result.pos = node->getDim();
                invariantBroken |= result.best == localX;
            }
            if (node->getLo() != null) {
                removeMaxLeaf(node->getLo(), node, pos, result);
            }
            if (node->getHi() != null) {
                removeMaxLeaf(node->getHi(), node, pos, result);
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
    const T& update(const Key& oldKey, const Key& newKey) {
        if (root == null) {
            return null;
        }
        const T& value = remove(oldKey);
        insert(newKey, value);
        return value;
    }

    /**
     * Get the number of key-value pairs in the tree.
     * @return the size
     */
  public:
    int size() {
        return size_;
    }

    /**
     * Removes all elements from the tree.
     */
  public:
    void clear() {
        size_ = 0;
        root = null;
        invariantBroken = false;
        modCount_++;
    }

    /**
     * Query the tree, returning all points in the axis-aligned rectangle between 'min' and 'max'.
     * @param min lower left corner of query
     * @param max upper right corner of query
     * @return all entries in the rectangle
     */
    public;
    KDIterator<Key, T> query(const Key& min, const Key& max) {
        return new KDIterator(this, min, max);
    }

  private:
    double distance(const Key& p1, const Key& p2) {
        return dist.dist(p1, p2);
    }

    /**
     * 1-nearest neighbor query.
     * @param center The point for which the nearest neighbors are requested
     * @return Nearest neighbor
     */
  public:
    KDEntryDist<Key, T> nnQuery(const Key& center) {
        if (root == nullptr) {
            return nullptr;
        }
        KDEntryDist<Key, T> candidate = new KDEntryDist(nullptr, Double.POSITIVE_INFINITY);
        rangeSearch1NN(root, center, candidate, Double.POSITIVE_INFINITY);
        return candidate;
    }

  private:
    double rangeSearch1NN(
        Node<Key, T>* node, const Key& center, KDEntryDist<Key, T> candidate, double maxRange) {
        int pos = node->getDim();
        if (node->getLo() != nullptr &&
            (center[pos] < node->getKey()[pos] || node->getHi() == nullptr)) {
            // go down
            maxRange = rangeSearch1NN(node->getLo(), center, candidate, maxRange);
            // refine result
            if (center[pos] + maxRange >= node->getKey()[pos]) {
                maxRange = addCandidate(node, center, candidate, maxRange);
                if (node->getHi() != nullptr) {
                    maxRange = rangeSearch1NN(node->getHi(), center, candidate, maxRange);
                }
            }
        } else if (node->getHi() != nullptr) {
            // go down
            maxRange = rangeSearch1NN(node->getHi(), center, candidate, maxRange);
            // refine result
            if (center[pos] <= node->getKey()[pos] + maxRange) {
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

  private:
    double addCandidate(
        Node<Key, T>* node,
        const Key& center,
        const KDEntryDist<Key, T>& candidate,
        double maxRange) {
        nDist1NN++;
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
    auto knnQuery(const Key& center, int k) {
        if (root == nullptr) {
            return std::vector<KDEntryDist<Key, T>>(0);
        }
        std::vector<KDEntryDist<Key, T>> candidates(k);
        rangeSearchKNN(root, center, candidates, k, std::numeric_limits<SCALAR>::infinity());
        return candidates;
    }

  private:
    double rangeSearchKNN(
        Node<Key, T>* node,
        const Key& center,
        std::vector<KDEntryDist<Key, T>>& candidates,
        int k,
        double maxRange) {
        int pos = node->getDim();
        if (node->getLo() != nullptr &&
            (center[pos] < node->getKey()[pos] || node->getHi() == nullptr)) {
            // go down
            maxRange = rangeSearchKNN(node->getLo(), center, candidates, k, maxRange);
            // refine result
            if (center[pos] + maxRange >= node->getKey()[pos]) {
                maxRange = addCandidate(node, center, candidates, k, maxRange);
                if (node->getHi() != nullptr) {
                    maxRange = rangeSearchKNN(node->getHi(), center, candidates, k, maxRange);
                }
            }
        } else if (node->getHi() != nullptr) {
            // go down
            maxRange = rangeSearchKNN(node->getHi(), center, candidates, k, maxRange);
            // refine result
            if (center[pos] <= node->getKey()[pos] + maxRange) {
                maxRange = addCandidate(node, center, candidates, k, maxRange);
                if (node->getLo() != nullptr) {
                    maxRange = rangeSearchKNN(node->getLo(), center, candidates, k, maxRange);
                }
            }
        } else {
            // leaf -> first (probably best) match!
            maxRange = addCandidate(node, center, candidates, k, maxRange);
        }
        return maxRange;
    }

  private:
    static const Comparator < KDEntryDist < ? >>
        compKnn = (KDEntryDist < ? > point1, KDEntryDist < ? > point2)->{
        double deltaDist = point1.dist() - point2.dist();
        return deltaDist < 0 ? -1 : (deltaDist > 0 ? 1 : 0);
    };

  private:
    double addCandidate(
        Node<Key, T>* node,
        const Key& center,
        std::vector<KDEntryDist<Key, T>>& candidates,
        int k,
        double maxRange) {
        nDistKNN++;
        // add ?
        double dist = distance(center, node->getKey());
        if (dist > maxRange) {
            // don't add if too far away
            return maxRange;
        }
        if (dist == maxRange && candidates.size() >= k) {
            // don't add if we already have enough equally good results.
            return maxRange;
        }
        KDEntryDist<Key, T> cand;
        if (candidates.size() >= k) {
            cand = candidates.remove(k - 1);
            cand.set(node, dist);
        } else {
            cand = new KDEntryDist(node, dist);
        }
        int insertionPos = Collections.binarySearch(candidates, cand, compKnn);
        insertionPos = insertionPos >= 0 ? insertionPos : -(insertionPos + 1);
        candidates.add(insertionPos, cand);
        return candidates.size() < k ? maxRange : candidates.get(candidates.size() - 1).dist();
    }

  private:
    static class KDQueryIteratorKNN<T> implements QueryIteratorKNN<PointEntryDist<T>> {
      private
        Iterator < ? extends PointEntryDist<T> > it;
      private
        final KDTree<T> tree;

      public
        KDQueryIteratorKNN(KDTree<T> tree, const Key& center, int k) {
            this.tree = tree;
            reset(center, k);
        }

      public
        boolean hasNext() {
            return it.hasNext();
        }

      public
        PointEntryDist<T> next() {
            return it.next();
        }

      public
        KDQueryIteratorKNN<T> reset(const Key& center, int k) {
            it = tree.knnQuery(center, k).iterator();
            return this;
        }
    }

    /**
     * Returns a printable list of the tree.
     * @return the tree as String
     */
    public : String
             toStringTree() {
        StringBuilder sb = new StringBuilder();
        if (root == null) {
            sb.append("empty tree");
        } else {
            toStringTree(sb, root, 0);
        }
        return sb.toString();
    }

  private:
    void toStringTree(StringBuilder sb, Node<Key, T>* node, int depth) {
        String prefix = "";
        for (int i = 0; i < depth; i++) {
            prefix += ".";
        }
        // sb.append(prefix + " d=" + depth + NL);
        prefix += " ";
        if (node->getLo() != null) {
            toStringTree(sb, node->getLo(), depth + 1);
        }
        sb.append(prefix + Arrays.toString(node->point()));
        sb.append(" v=" + node->value());
        sb.append(" l/r=");
        sb.append(node->getLo() == null ? null : Arrays.toString(node->getLo().point()));
        sb.append("/");
        sb.append(node->getHi() == null ? null : Arrays.toString(node->getHi().point()));
        sb.append(NL);
        if (node->getHi() != null) {
            toStringTree(sb, node->getHi(), depth + 1);
        }
    }

  public:
    String toString() {
        return "KDTree;size=" + size_ + ";DEBUG=" + DEBUG +
            ";DistFn=" + PointDistanceFunction.getName(dist) +
            ";center=" + (root == null ? "null" : Arrays.toString(root.getKey()));
    }

  public:
    KDStats getStats() {
        KDStats s = new KDStats(this);
        if (root != null) {
            root.checkNode(s, 0);
        }
        return s;
    }

    /**
     * Statistics container class.
     */
  public:
    static class KDStats extends Stats {
      public KDStats(KDTree<?> tree) {
            super(tree.nDist1NN + tree.nDistKNN, tree.nDist1NN, tree.nDistKNN);
        }
    }

    public : int
             getDims() {
        return dims_;
    }

  public:
    QueryIterator<PointEntry<T>> iterator() {
        if (root == null) {
            return query(new double[dims_], new double[dims_]);
        }
        // return query(root.);
        // TODO
        throw new UnsupportedOperationException();
    }

  public:
    KDEntryDist<T> query1NN(const Key& center) {
        return nnQuery(center);
    }

  public:
    QueryIteratorKNN<PointEntryDist<T>> queryKNN(const Key& center, int k) {
        return new KDQueryIteratorKNN<>(this, center, k);
    }

  public:
    int getNodeCount() {
        return getStats().getNodeCount();
    }

  public:
    int getDepth() {
        return getStats().getMaxDepth();
    }

    Node<Key, T>* getRoot() {
        return root;
    }
};

}  // namespace tinspin

#endif