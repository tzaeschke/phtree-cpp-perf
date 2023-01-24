// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef TINSPIN_QUADTREE_POINT_H
#define TINSPIN_QUADTREE_POINT_H

#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include <iostream>
#include <cmath>
#include <stack>
#include <unordered_set>

namespace tinspin {

using namespace improbable::phtree;

/**
 * A simple MX-quadtree implementation with configurable maximum depth, maximum nodes size, and
 * (if desired) automatic guessing of root rectangle.
 *
 * This version of the quadtree stores for each node only the center point and the
 * distance (radius) to the edges.
 * This reduces space requirements but increases problems with numerical precision.
 * Overall it is more space efficient and slightly faster.
 *
 * @author ztilmann
 *
 * @param <T> Value type.
 */

namespace {

static const bool DEBUG = false;

template <typename KEy, typename T>
class QREntry;


static const double EPS_MUL = 1.000000001;

template<typename Key>
static bool isPointEnclosed(const Key& point,
                               const Key& min, const Key& max) {
    for (size_t d = 0; d < min.size(); ++d) {
        if (point[d] < min[d] || point[d] > max[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
static bool isPointEnclosed(const Key& point,
                               const Key& center, double radius) {
    for (size_t d = 0; d < center.size(); ++d) {
        if (point[d] < center[d]-radius || point[d] > center[d]+radius) {
            return false;
        }
    }
    return true;
}

template<typename Key>
static bool isPointEqual(const Key& p1, const Key& p2) {
    for (size_t d = 0; d < p1.size(); ++d) {
        if (p1[d] != p2[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
static bool isRectEqual(const Key& p1L, const Key& p1U, const Key& p2L, const Key& p2U) {
    return isPointEqual(p1L, p2L) && isPointEqual(p1U, p2U);
}

template<typename Key, typename T>
static bool isRectEqual(QREntry<Key, T> e, const Key& keyL, const Key& keyU) {
    return isRectEqual(e.lower(), e.upper(), keyL, keyU);
}

template<typename Key>
static bool overlap(const Key& min, const Key& max, const Key& min2, const Key& max2) {
    for (size_t d = 0; d < min.size(); ++d) {
        if (max[d] < min2[d] || min[d] > max2[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool overlap(const Key& min, const Key& max, const Key& center, double radius) {
    for (size_t d = 0; d < min.size(); ++d) {
        if (max[d] < center[d]-radius || min[d] > center[d]+radius) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool overlap(const Key& center, double radius,
                       const Key& center2, double radius2) {
    for (size_t d = 0; d < center.size(); ++d) {
        if (center[d]+radius < center2[d]-radius2 || center[d]-radius > center2[d]+radius2) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool isRectEnclosed(const Key& minEnclosed, const Key& maxEnclosed,
                              const Key& minOuter, const Key& maxOuter) {
    for (size_t d = 0; d < minOuter.size(); ++d) {
        if (maxOuter[d] < maxEnclosed[d] || minOuter[d] > minEnclosed[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool isRectEnclosed(const Key& minEnclosed, const Key& maxEnclosed,
                              const Key& centerOuter, double radiusOuter) {
    for (size_t d = 0; d < centerOuter.size(); ++d) {
        double radOuter = radiusOuter;
        if ((centerOuter[d]+radOuter) < maxEnclosed[d] ||
            (centerOuter[d]-radOuter) > minEnclosed[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool isRectEnclosed(const Key& centerEnclosed, double radiusEnclosed,
                              const Key& centerOuter, double radiusOuter) {
    for (size_t d = 0; d < centerOuter.size(); ++d) {
        double radOuter = radiusOuter;
        double radEncl = radiusEnclosed;
        if ((centerOuter[d]+radOuter) < (centerEnclosed[d]+radEncl) ||
            (centerOuter[d]-radOuter) > (centerEnclosed[d]-radEncl)) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static double distance(const Key& p1, const Key& p2) {
    double dist = 0;
    for (size_t i = 0; i < p1.size(); i++) {
        double d = p1[i]-p2[i];
        dist += d * d;
    }
    return std::sqrt(dist);
}

/**
	 * Calculates distance to center point of rectangle.
	 * @param p point
	 * @param rMin rectangle min
	 * @param rMax rectangle max
	 * @return distance to center point
 */
template<typename Key>
 static double distToRectCenter(const Key& p, const Key& rMin, const Key& rMax) {
    double dist = 0;
    for (size_t i = 0; i < p.size(); i++) {
        double d = (rMin[i]+rMax[i])/2. - p[i];
        dist += d * d;
    }
    return std::sqrt(dist);
}

/**
	 * Calculates distance to center point of rectangle.
	 * @param p point
	 * @param e rectangle
	 * @return distance to center point
 */
template<typename Key, typename T>
 static double distToRectCenter(const Key& p, QREntry<Key, T> e) {
    return distToRectCenter(p, e.lower(), e.upper());
}

/**
	 * Calculates distance to the edge of rectangle.
	 * @param p point
	 * @param rMin rectangle min
	 * @param rMax rectangle max
	 * @return distance to edge
 */
template<typename Key>
 double distToRectEdge(const Key& center, const Key& rLower, const Key& rUpper) {
    double dist = 0;
    for (size_t i = 0; i < center.size(); i++) {
        double d = 0;
        if (center[i] > rUpper[i]) {
            d = center[i] - rUpper[i];
        } else if (center[i] < rLower[i]) {
            d = rLower[i] - center[i];
        }
        dist += d*d;
    }
    return std::sqrt(dist);
}

/**
	 * Calculates distance to edge of rectangle.
	 * @param p point
	 * @param e rectangle
	 * @return distance to edge point
 */
template<typename Key, typename T>
 static double distToRectEdge(const Key& p, QREntry<Key, T> e) {
    return distToRectEdge(p, e.lower(), e.upper());
}

/**
	 * Calculates distance to the edge of a node.
	 * @param point the point
	 * @param nodeCenter the center of the node
	 * @param nodeRadius radius of the node
	 * @return distance to edge of the node or 0 if the point is inside the node
 */
template<typename Key>
static double distToRectNode(const Key& point, const Key& nodeCenter, double nodeRadius) {
    double dist = 0;
    for (size_t i = 0; i < point.size(); i++) {
        double d = 0;
        if (point[i] > nodeCenter[i]+nodeRadius) {
            d = point[i] - (nodeCenter[i]+nodeRadius);
        } else if (point[i] < nodeCenter[i]-nodeRadius) {
            d = nodeCenter[i]-nodeRadius - point[i];
        }
        dist += d*d;
    }
    return std::sqrt(dist);
}

/**
 * Statistics container class.
 */
class QStats {
  public:
    int dims;
    int nEntries = 0;
    int nNodes = 0;
    int minLevel = std::numeric_limits<int>::max();
    int maxLevel = -1;
    int maxDepth = 0;
    double sumLevel;
    int maxNodeSize = -1;
    int nLeaf;
    int nInner;
    long nDistCalc;
    long nDistCalc1NN;
    long nDistCalcKNN;

    QStats(long nDistCalc = 0, long nDistCalc1NN = 0, long nDistCalcKNN = 0) {
        this->nDistCalc = nDistCalc;
        this->nDistCalc1NN = nDistCalc1NN;
        this->nDistCalcKNN = nDistCalcKNN;
    }

    void toString() {
        std::cout << "dims=" << dims << ";nEntries=" << nEntries << ";nNodes=" << nNodes
                  << ";nLeaf=" << nLeaf << ";nInner=" << nInner << ";minLevel=" << minLevel
                  << ";maxLevel=" << maxLevel << ";avgLevel=" << (sumLevel / nEntries)
                  << ";maxNodeSize=" << maxNodeSize;
    }
};

template <typename Key, typename T>
class QEntry {
    const Key& point_;
    T value_;

  public:
    QEntry(const Key& key, const T& value) : point_{key}, value_{value} {}

    const Key& point() const {
        return point_;
    }

    const T& value() const {
        return value_;
    }

    bool enclosedBy(const Key& min, const Key& max) const {
        return isPointEnclosed(point_, min, max);
    }

    bool enclosedBy(const Key& center, double radius) const {
        return isPointEnclosed(point_, center, radius);
    }

    bool isExact(const QEntry<Key, T>& e) const {
        return isPointEqual(point_, e.point());
    }

    void setKey(const Key& newPoint) {
        point_ = newPoint;
    }
};

template <typename Key, typename T>
class QEntryDist : public QEntry<Key, T> {  // TODO inheritance???
    double distance_;

  public:
    QEntryDist(const QEntry<Key, T>& e, double dist) : QEntry<Key, T>(e.point(), e.value()) {
        distance_ = dist;
    }

    double dist() const {
        return distance_;
    }
};

/**
 * Node class for the quadtree.
 *
 * @author ztilmann
 *
 * @param <T> Value type.
 */
template<typename Key, typename T>
class QNode {
    Key center_;
    double radius_;
    // null indicates that we have sub-nopde i.o. values
    std::vector<QEntry<Key, T>> values_;
    std::vector<QNode<Key, T>*> subs_;
    bool is_leaf_;

  public:
    QNode(const Key& center, double radius) : center_{center}, radius_{radius}, values_{}, subs_{}, is_leaf_{true} {
        values_.reserve(2);
    }

    QNode(const Key& center, double radius, QNode<Key, T>* subNode) : center_{center}, radius_{radius}, values_{}, subs_{}, is_leaf_{false} {
        subs_.emplace_back(subNode);
    }

    QNode<Key, T>* tryPut(const QEntry<Key, T>& e, size_t maxNodeSize, bool enforceLeaf) {
        if (DEBUG && !e.enclosedBy(center_, radius_)) {
            std::cerr << "ERROR: e=" << e.point() << " center/radius=" << center_ << "/" << radius_ << std::endl;
        }

        // traverse subs?
        if (!isLeaf()) {
            return getOrCreateSub(e);
        }

        // add if:
        // a) we have space
        // b) we have maxDepth
        // c) elements are equal (work only for n=1, avoids splitting
        //    in cases where splitting won't help. For n>1 the
        //    local limit is (temporarily) violated.
        if (values_.size() < maxNodeSize || enforceLeaf || e.isExact(values_[0])) {
            values_.emplace_back(e);
            return nullptr;
        }

        // split
        std::vector<QEntry<Key, T>> vals = std::move(values_);  // TODO avoid move and erase later?
        values_.clear(); // = nullptr;
        values_.shrink_to_fit();
        assert(subs_.empty());
        for (size_t i = 0; i < vals.size(); i++) {
            QEntry<Key, T>& e2 = vals[i];
            QNode<Key, T>* sub = getOrCreateSub(e2);
            while (sub != nullptr) {
                // This may recurse if all entries fall
                // into the same subnode
                sub = sub->tryPut(e2, maxNodeSize, false);
            }
        }
        return getOrCreateSub(e);
    }

  private:
    QNode<Key, T>* getOrCreateSub(const QEntry<Key, T>& e) {
        QNode<Key, T>* n = findSubNode(e.point());
        if (n == nullptr) {
            n = createSubForEntry(e.point());
            subs_.emplace_back(n);
        }
        return n;
    }

    QNode<Key, T>* createSubForEntry(const Key& p) const {
        Key centerSub{};
        // This ensures that the subsnodes completely cover the area of
        // the parent node.
        double radiusSub = radius_ / 2.0;
        for (size_t d = 0; d < center_.size(); ++d) {
            if (p[d] >= center_[d]) {
                centerSub[d] = center_[d] + radiusSub;
            } else {
                centerSub[d] = center_[d] - radiusSub;
            }
        }
        return new QNode<Key, T>(std::move(centerSub), radiusSub);
    }

    /**
	 * The subnode position has reverse ordering of the point's
	 * dimension ordering. Dimension 0 of a point is the highest
	 * ordered bit in the position.
	 * @param p point
	 * @return subnode position
     */
    QNode<Key, T>* findSubNode(const Key& p) const {
        for (size_t i = 0; i < subs_.size(); i++) {
            QNode<Key, T>* n = subs_[i];
            if (isPointEnclosed(p, n->center_, n->radius_)) {
                return n;
            }
        }
        return nullptr;
    }

  public:
    size_t remove(QNode<Key, T>* parent, const Key& key, int maxNodeSize) {
        if (!is_leaf_) {
            QNode<Key, T>* sub = findSubNode(key);
            if (sub != nullptr) {
                return sub->remove(this, key, maxNodeSize);
            }
            return 0;
        }

        size_t n = 0;
        for (size_t i = 0; i < values_.size(); i++) {
            QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.point(), key)) {
                values_.erase(values_.begin + i);
                ++n;
            }
        }
        // TODO provide threshold for re-insert
        // i.e. do not always merge.
        if (n > 0 && parent != nullptr) {
            parent->checkAndMergeLeafNodes(maxNodeSize);
        }
        return n;
    }

    size_t remove(QNode<Key, T>* parent, const Key& key, int maxNodeSize, const T& value) {
        if (values_ == nullptr) {
            QNode<Key, T>* sub = findSubNode(key);
            if (sub != nullptr) {
                return sub->remove(this, key, maxNodeSize);
            }
            return 0;
        }

        size_t n = 0;
        for (size_t i = 0; i < values_.size(); i++) {
            QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.point(), key) && e.value() == value) {
                values_.remove(i);
                ++n;
            }
        }
        // TODO provide threshold for re-insert
        // i.e. do not always merge.
        if (n > 0 && parent != nullptr) {
            parent->checkAndMergeLeafNodes(maxNodeSize);
        }
        return n;
    }

     QEntry<Key, T>* update(
        QNode<Key, T>* parent,
        const Key& keyOld,
        const Key& keyNew,
        int maxNodeSize,
        bool& requiresReinsert,
        int currentDepth,
        int maxDepth) {
        if (values_ == nullptr) {
            QNode<Key, T>* sub = findSubNode(keyOld);
            if (sub == nullptr) {
                return nullptr;
            }
            QEntry<Key, T>& ret = sub->update(
                this, keyOld, keyNew, maxNodeSize, requiresReinsert, currentDepth + 1, maxDepth);
            if (ret != nullptr && requiresReinsert &&
                isPointEnclosed(ret.point(), center_, radius_)) {
                requiresReinsert = false;
                auto* r = this;  // TODO backport, r is always a QNode
                while (r != nullptr) {
                    r = r->tryPut(ret, maxNodeSize, currentDepth++ > maxDepth);
                }
            }
            return ret;
        }

        for (size_t i = 0; i < values_.size(); i++) {
            QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.point(), keyOld)) {
                values_.remove(i);
                e.setKey(keyNew);
                if (isPointEnclosed(keyNew, center_, radius_)) {
                    // reinsert locally;
                    values_.add(e);
                    requiresReinsert = false;
                } else {
                    requiresReinsert = true;
                    // TODO provide threshold for re-insert
                    // i.e. do not always merge.
                    if (parent != nullptr) {
                        parent->checkAndMergeLeafNodes(maxNodeSize);
                    }
                }
                return e;
            }
        }
        requiresReinsert = false;
        return nullptr;
    }

  private:
    void checkAndMergeLeafNodes(int maxNodeSize) {
        // check
        int nTotal = 0;
        for (size_t i = 0; i < subs_.size(); i++) {
            if (!subs_[i]->is_leaf_) {
                // can't merge directory nodes.
                return;
            }
            nTotal += subs_[i]->values_.size();
            if (nTotal > maxNodeSize) {
                // too many children
                return;
            }
        }

        // okay, let's merge
        assert(values_.empty());
        values_.reserve(nTotal);
        //values_ = new ArrayList<>(nTotal);
        for (size_t i = 0; i < subs_.size(); i++) {
            values_.insert(subs_[i]->values_.begin(), subs_[i]->values_.end());
        }
        //subs = nullptr;
        subs_.clear();
        subs_.shrink_to_fit();
    }

  public:
    const Key& getCenter() {
        return center_;
    }

    double getRadius() {
        return radius_;
    }

    const QEntry<Key, T>* getExact(const Key& key) const {
        if (!is_leaf_) {
            QNode<Key, T>* sub = findSubNode(key);
            if (sub != nullptr) {
                return sub->getExact(key);
            }
            return nullptr;
        }

        for (size_t i = 0; i < values_.size(); ++i) {
            const QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.point(), key)) {
                return &e;
            }
        }
        return nullptr;
    }

    const QEntry<Key, T>* getExact(const Key& key, const T& value) const {
        if (!is_leaf_) {
            QNode<Key, T>* sub = findSubNode(key);
            if (sub != nullptr) {
                return sub->getExact(key);
            }
            return nullptr;
        }

        for (size_t i = 0; i < values_.size(); ++i) {
            const QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.point(), key) && e.value() == value) {
                return &e;
            }
        }
        return nullptr;
    }

    auto& getEntries() {
        return values_;
    }

//    auto getChildIterator() {
//        if (values_ != nullptr) {
//            return values_.iterator();
//        }
//        return subs_.iterator();
//    }
    auto getChildNodeIterator() const {
        assert(!subs_.empty());
        assert(!isLeaf());
        return subs_.begin();
    }

    auto getChildEntryIterator() {
        assert(isLeaf());
        return values_.begin();
    }

//    public: String toString() {
//        return "center/radius=" + Arrays.toString(center) + "/" + radius + " " +
//            System.identityHashCode(this);
//    }

    void checkNode(QStats s, QNode<Key, T>* parent, int depth) {
        if (depth > s.maxDepth) {
            s.maxDepth = depth;
        }
        s.nNodes++;

        if (parent != nullptr) {
            if (!isRectEnclosed(
                    center_, radius_, parent->center_, parent->radius * EPS_MUL)) {
                for (size_t d = 0; d < center_.size(); ++d) {
                    //					if ((centerOuter[d]+radiusOuter) / (centerEnclosed[d]+radiusEnclosed) < 0.9999999 || 							(centerOuter[d]-radiusOuter) / (centerEnclosed[d]-radiusEnclosed) > 1.0000001) { 						return false;
                    //					}
                    std::cout <<
                        "Outer: " <<  parent->radius_ <<  " " <<  parent->center_ << std::endl;
                    std::cout << "Child: " <<  radius_ <<  " " <<  center_ << std::endl;
                    std::cout <<
                        (parent->center_[d] +  parent->radius_) <<  " vs " <<  (center_[d] +  radius_) << std::endl;
                    std::cout <<
                        "r=" <<  (parent->center_[d] +  parent->radius_) / (center_[d] +  radius_) << std::endl;
                    std::cout <<
                        (parent->center_[d] - parent->radius_) <<  " vs " <<  (center_[d] - radius_) << std::endl;
                    std::cout <<
                        "r=" <<  (parent->center_[d] - parent->radius_) / (center_[d] - radius_) << std::endl;
                }
                assert(false);
            }
        }
        if (values_ != nullptr) {
            for (size_t i = 0; i < values_.size(); i++) {
                QEntry<Key, T>& e = values_[i];
                if (!isPointEnclosed(e.point(), center_, radius_ * EPS_MUL)) {
                    std::cout << "Node: " <<  radius_ <<  " " <<  center_ << std::endl;
                    std::cout << "Child: " <<  e.point() << std::endl;
                    for (size_t d = 0; d < center_.size(); ++d) {
                        //						if ((centerOuter[d]+radiusOuter) / (centerEnclosed[d]+radiusEnclosed) < 0.9999999 ||
                        //								(centerOuter[d]-radiusOuter) / (centerEnclosed[d]-radiusEnclosed) > 1.0000001) { 							return false;
                        //						}
                        std::cout << "min/max for " <<  d << std::endl;
                        std::cout <<
                            "min: " <<  (center_[d] - radius_) <<  " vs " <<  (e.point()[d]) << std::endl;
                        std::cout << "r=" <<  (center_[d] - radius_) / (e.point()[d]) << std::endl;
                        std::cout <<
                            "max: " <<  (center_[d] <<  radius_) <<  " vs " <<  (e.point()[d]) << std::endl;
                        std::cout << "r=" <<  (center_[d] <<  radius_) / (e.point()[d]) << std::endl;
                    }
                    assert(false);
                }
            }
            if (subs_ != nullptr) {
                assert(false);
            }
        } else {
            for (size_t i = 0; i < subs_.size(); i++) {
                QNode<Key, T>& n = subs_[i];
                n.checkNode(s, this, depth + 1);
            }
        }
    }

    bool isLeaf() const noexcept {
        return is_leaf_;
    }

    auto& getChildNodes() {
        return subs_;
    }
};

template <typename Key, typename T>
class QIteratorBase {
  public:
    explicit QIteratorBase() noexcept : node_{nullptr} {}

    inline auto& operator*() const noexcept {
        assert(node_ != nullptr);
        return entry_->value();
    }

    inline auto* operator->() const noexcept {
        assert(node_ != nullptr);
        return &entry_->value();
    }

    inline friend bool operator==(
        const QIteratorBase<Key, T>& left, const QIteratorBase<Key, T>& right) noexcept {
        // TODO compare stack status left/right/key
        // TODO do not compare node, compare only Entry*...!!!!!!!!!!!!!!!
        return left.node_ == right.node_ && left.entry_ == right.entry_;
    }

    inline friend bool operator!=(
        const QIteratorBase<Key, T>& left, const QIteratorBase<Key, T>& right) noexcept {
        return left.node_ != right.node_ || left.entry_ == right.entry_;
    }

        auto _node()  const noexcept {
            return node_;
        }

                auto _entry()  const noexcept {
                    return entry_;
                }

  protected:
    bool IsEnd() const noexcept {
        return this->_node() == nullptr;
    }

    void SetFinished() noexcept {
        node_ = nullptr;
    }

        void SetCurrentNode(QNode<Key, T>* node) noexcept {
            node_ = node;
        }

                void SetCurrentResult(QEntry<Key, T>* entry) noexcept {
            entry_ = entry;
                }

  protected:
    QNode<Key, T>* node_ = nullptr;
    QEntry<Key, T>* entry_ = nullptr;
};


template <typename Key, typename T, typename FILTER>
class QIterator : public QIteratorBase<Key, T> {
    using IterNodeT = decltype(std::vector<QNode<Key, T>*>().begin());
    using IterEntryT = decltype(std::vector<QEntry<Key, T>>().begin());
    using EntryInnerT = std::pair<QNode<Key, T>*, IterNodeT>;
    const QNode<Key, T>* root_;
    std::stack <EntryInnerT> stack_; // TODO backport, Why is this a Deque????????
    IterEntryT iter_leaf_;
    Key min;
    Key max;
    FILTER filter_;

  public:
    template <typename F = FilterNoOp>
    QIterator(QNode<Key, T>* root, const Key& min, const Key& max, F&& filter = F())
    : QIteratorBase<Key, T>()
    , root_{root}
    , stack_{}
    , min(min)
    , max(max)
    , filter_{std::forward<F>(filter)} {
            if (root != nullptr) {
            this->SetCurrentNode(root);
            if (root->isLeaf()) {
                iter_leaf_ = root->getChildEntryIterator();
            } else {
                stack_.emplace(std::make_pair(root, root->getChildNodeIterator()));
            }
            findNext();
            } else {
            this->SetFinished();
            }
    }

    QIterator& operator++() noexcept {
        assert(!this->IsEnd());
        findNext();
        return *this;
    }

    QIterator operator++(int) noexcept {
        assert(!this->IsEnd());
        QIterator iterator(*this);
        ++(*this);
        return iterator;
    }

  private: void findNext() {
        assert(!this->IsEnd());
        do {
            // check current leaf
            auto * current = this->node_;
            if (current->isLeaf()) {
                while (iter_leaf_ != current->getEntries().end()) {
                    QEntry<Key, T>& e = *iter_leaf_;
                    ++iter_leaf_;
                    if (e.enclosedBy(min, max) && filter_.IsEntryValid(e.point(), e.value())) {
                        this->entry_ = &e;
                        return;
                    }
                }
            }
            if (stack_.empty()) {
                break;
            }
            // traverse inner nodes
            auto& ee = stack_.top();
            auto &it = ee.second;
            while (it != ee.first->getChildNodes().end()) {
                auto* node = *it;
                this->node_ = node;
                if (overlap(min, max, node->getCenter(), node->getRadius())) {
                    if (node->isLeaf()) {
                        iter_leaf_ = node->getChildEntryIterator();
                    } else {
                        auto it2 = node->getChildNodes().begin();
                        stack_.push(it2);
                    }
                    break;
                }
                ++it;
            }
            stack_.pop();
        } while (!stack_.empty());
        this->SetFinished();
    }
};

template<typename Key, typename T>
class QIteratorFind : public QIteratorBase<Key, T> {
    using IterNodeT = decltype(std::vector<QNode<Key, T>>().begin());
    using IterEntryT = decltype(std::vector<QEntry<Key, T>>().begin());
    using EntryInnerT = std::pair<QNode<Key, T>*, IterNodeT>;

  public:
    template <typename F = FilterNoOp>
    QIteratorFind(QNode<Key, T>* root, QEntry<Key, T>* entry) : QIteratorBase<Key, T>() {
        this->SetCurrentNode(root);
        this->entry_ = entry;
    }

    QIteratorFind& operator++() noexcept {
        assert(!this->IsEnd());
        this->SetFinished(); // TODO ....
        return *this;
    }

    QIteratorFind operator++(int) noexcept {
        assert(!this->IsEnd());
        QIterator iterator(*this);
        ++(*this);
        return iterator;
    }
};

template <typename Key, typename T>
class QIteratorKnn : public QIteratorBase<Key, T> {
    using Candidates = std::vector<QEntryDist<Key, T>>;
    using CandidatesIter = decltype(Candidates{}.begin());
  public:
    QIteratorKnn(Candidates&& result) noexcept
    : QIteratorBase<Key, T>(), result_{std::move(result)}, iter_{result_.begin()} {
        if (iter_ != result_.end()) {
            this->SetCurrentResult(iter_->_node());
        } else {
            this->SetFinished();
        }
    }

    QIteratorKnn& operator++() noexcept {
        assert(!this->IsEnd());
        findNext();
        return *this;
    }

    QIteratorKnn operator++(int) noexcept {
        QIteratorKnn iterator(*this);
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

}

template<typename T>
class QuadTree {
    static const int MAX_DEPTH = 50;
    using Key = PhPointD<3>;
    using QueryBox = PhBox<Key{}.size(), double>;

    static const int DEFAULT_MAX_NODE_SIZE = 10;

    const size_t dims;
    const size_t maxNodeSize;
    QNode<Key, T>* root_ = nullptr;
    size_t size_ = 0;

  public:
    using KeyInternal = Key;

    QuadTree(size_t dims = 3, size_t maxNodeSize = DEFAULT_MAX_NODE_SIZE) : dims{dims}, maxNodeSize{maxNodeSize} {
        if (DEBUG) {
            std::cout << "Warning: DEBUG enabled" << std::endl; // TODO
        }
    }


//    QuadTree<T> create(size_t dims, size_t maxNodeSize, const Key& center, double radius) {
//        QuadTree<T> t = new QuadTree<>(dims, maxNodeSize);
//        if (radius <= 0) {
//            throw new IllegalArgumentException("Radius must be > 0 but was " + radius);
//        }
//        t.root = new QNode<Key, T>(Arrays.copyOf(center, center.length), radius);
//        return t;
//    }

    /**
	 * Insert a key-value pair.
	 * @param key the key
	 * @param value the value
     */
    template <typename T2>
    void emplace(const Key& key, T2&& value) {
        size_++;
        QEntry<Key, T> e(key, std::forward<T2>(value)); // TODO std::move into node
        if (root_ == nullptr) {
            initializeRoot(key);
        }
        ensureCoverage(e);
        auto* r = root_;
        int depth = 0;
        while (r != nullptr) {  // TODO backport, r is always a QNode
            r = r->tryPut(e, maxNodeSize, depth++ > MAX_DEPTH);
        }
    }

    void insert(const Key& key, const T& value) {
        size_++;
        QEntry<Key, T> e(key, value); // TODO std::move into node
        if (root_ == nullptr) {
            initializeRoot(key);
        }
        ensureCoverage(e);
        auto* r = root_;
        int depth = 0;
        while (r != nullptr) {  // TODO backport, r is always a QNode
            r = r->tryPut(e, maxNodeSize, depth++ > MAX_DEPTH);
        }
    }

  private:
    void initializeRoot(const Key& key) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (size_t d = 0; d < dims; d++) {
            lo = lo > key[d] ? key[d] : lo;
            hi = hi < key[d] ? key[d] : hi;
        }
        if (lo == 0 && hi == 0) {
            hi = 1.0;
        }
        double maxDistOrigin = std::abs(hi) > std::abs(lo) ? hi : lo;
        maxDistOrigin = std::abs(maxDistOrigin);
        // no we use (0,0)/(+-maxDistOrigin*2,+-maxDistOrigin*2) as root.

        // HACK: To avoid precision problems, we ensure that at least the initial
        // point is not exactly on the border of the quadrants:
        maxDistOrigin *= EPS_MUL * EPS_MUL;
        Key center{};
        for (size_t d = 0; d < dims; d++) {
            center[d] = key[d] > 0 ? maxDistOrigin : -maxDistOrigin;
            //			max[d] = key[d] < 0 ? 0 : (maxDistOrigin*2);
        }
        root_ = new QNode<Key, T>(center, maxDistOrigin);
    }

  public:
    /**
	 * Check whether a given key exists.
	 * @param key the key to check
	 * @return true iff the key exists
     */
    size_t count(const Key& key) const {
        if (root_ == nullptr) {
            return 0;
        }
        return root_->getExact(key) != nullptr; // TODO
    }

    /**
	 * Get the value associates with the key.
	 * @param key the key to look up
	 * @return the value for the key or 'nullptr' if the key was not found
     */
    auto find(const Key& key) const {
        if (root_ == nullptr) {
            return QIteratorFind<Key, T>(nullptr, nullptr);
        }
        QEntry<Key, T>* e = const_cast<QEntry<Key, T>*>(root_->getExact(key));
//        return e == nullptr ? nullptr : e.value();

        return QIteratorFind<Key, T>(root_, e);
    }

    auto find(const Key& key, const T& value) const {
        if (root_ == nullptr) {
            return QIteratorFind<Key, T>(nullptr, nullptr);
        }
        QEntry<Key, T>* e = const_cast<QEntry<Key, T>*>(root_->getExact(key, value));
//        return e == nullptr ? nullptr : e.value();
        // TODO avoid using root
        return QIteratorFind<Key, T>(root_, e);
    }

    /**
	 * Remove a key.
	 * @param key key to remove
	 * @return the value associated with the key or 'nullptr' if the key was not found
     */
  size_t erase(const Key& key) {
        if (root_ == nullptr) {
            if (DEBUG) {
                std::cerr <<"Failed remove 1: " << key << std::endl;
            }
            return 0;
        }
        size_t n = root_->remove(nullptr, key, maxNodeSize);
        if (n == 0) {
            if (DEBUG) {
                std::cerr << "Failed remove 2: " << key << std::endl;
            }
            return 0;
        }
        size_-= n;
        return n;
    }

  size_t erase(const Key& key, const T& value) {
        if (root_ == nullptr) {
            if (DEBUG) {
                std::cerr <<"Failed remove 1: " << key << std::endl;
            }
            return 0;
        }
        size_t n = root_->remove(nullptr, key, maxNodeSize);
        if (n == 0) {
            if (DEBUG) {
                std::cerr << "Failed remove 2: " << key << std::endl;
            }
            return 0;
        }
        size_-= n;
        return n;
    }

    /**
	 * Reinsert the key.
	 * @param oldKey old key
	 * @param newKey new key
	 * @return the value associated with the key or 'nullptr' if the key was not found.
     */
  size_t relocate(const Key& oldKey, const Key& newKey) {
        if (root_ == nullptr) {
            return 0;
        }
        bool requiresReinsert = false;
        QEntry<Key, T>* e =
            root_->update(nullptr, oldKey, newKey, maxNodeSize, requiresReinsert, 0, MAX_DEPTH);
        if (e == nullptr) {
            // not found
            if (DEBUG) {
                std::cout << "Failed reinsert 1: " << oldKey << "/" << newKey << std::endl;
            }
            return 0;
        }
        if (requiresReinsert) {
            if (DEBUG) {
                std::cout << "Failed reinsert 2: " << oldKey << "/" << newKey << std::endl;
            }
            // does not fit in root node...
            ensureCoverage(e);
            auto* r = root_;
            int depth = 0;
            while (r != nullptr) {  // TODO backport, r is always a QNode
                r = r->tryPut(e, maxNodeSize, depth++ > MAX_DEPTH);
            }
        }
        return 1;
    }

    /**
	 * Ensure that the tree covers the entry.
	 * @param e Entry to cover.
     */
    private: void ensureCoverage(const QEntry<Key, T>& e) {
        const Key& p = e.point();
        while (!e.enclosedBy(root_->getCenter(), root_->getRadius())) {
            const Key& center = root_->getCenter();
            double radius = root_->getRadius();
            Key center2{};
            double radius2 = radius * 2;
            // TODO use DIM?
            for (size_t d = 0; d < center.size(); d++) {
                if (p[d] < center[d] - radius) {
                    center2[d] = center[d] - radius;
                    // root will end up in upper quadrant in this
                    // dimension
                } else {
                    // extend upwards, even if extension unnecessary for this dimension.
                    center2[d] = center[d] + radius;
                }
            }
            if (DEBUG && !isRectEnclosed(center, radius, center2, radius2)) {
                std::cout <<
                    "e=" << e.point() <<
                    " center/radius=" << center2 << "/" << radius << std::endl;
            }
            root_ = new QNode<Key, T>(center2, radius2, root_);
        }
    }

    /**
	 * Get the number of key-value pairs in the tree.
	 * @return the size
     */
    public: int size() {
        return size_;
    }

    /**
	 * Removes all elements from the tree.
     */
    void clear() {
        size_ = 0;
        root_ = nullptr;
    }

    /**
	 * Query the tree, returning all points in the axis-aligned rectangle between 'min' and 'max'.
	 * @param min lower left corner of query
	 * @param max upper right corner of query
	 * @return all entries in the rectangle
     */
//    public: auto query(const Key& min, const Key& max) {
//        return QIterator<Key, T>(root_, min, max);
//    }

    template <typename DISTANCE, typename FILTER = FilterNoOp>
    std::vector<QEntryDist<Key, T>>
    knnQuery(const Key& center, size_t k,
             DISTANCE&& distance_function = DISTANCE(),
             FILTER&& filter = FILTER()) const {
        if (root_ == nullptr) {
            return std::vector<QEntryDist<Key, T>>{};
        }
        auto comp = [&center](const QEntry<Key, T>& point1, const QEntry<Key, T>& point2) {
//            double deltaDist =
//                distance(center, point1.point()) - distance(center, point2.point());
//            return deltaDist < 0 ? -1 : (deltaDist > 0 ? 1 : 0);
            return
                distance(center, point1.point()) < distance(center, point2.point());
        };
        double distEstimate = distanceEstimate(root_, center, k, comp);
        std::vector<QEntryDist<Key, T>> candidates{};
        candidates.reserve(k);
        while (candidates.size() < k) {
            candidates.clear();
            rangeSearchKNN(root_, center, candidates, k, distEstimate);
            distEstimate *= 2;
        }
        return candidates;
    }

    private:
      template <typename COMP>
      double distanceEstimate(QNode<Key, T>* node, const Key& point, size_t k, const COMP& comp) {
        if (node->isLeaf()) {
            // This is a leaf that would contain the point.
            size_t n = node->getEntries().size();
            // Create a copy!
            std::vector<QEntry<Key, T>> data(node->getEntries()); // TODO this is bad!!!! -> backport?
            std::sort(data.begin(), data.end(), comp);
            size_t pos = n < k ? n : k;
            double dist = distance(point, data[pos - 1].point());
            if (n < k) {
                // scale search dist with dimensions.
                dist = dist * std::pow(k / (double)n, 1 / (double)dims);
            }
            if (dist <= 0.0) {
                return node->getRadius();
            }
            return dist;
        } else {
            auto& nodes = node->getChildNodes();
            for (size_t i = 0; i < nodes.size(); i++) {
                QNode<Key, T>* sub = nodes[i];
                if (isPointEnclosed(point, sub->getCenter(), sub->getRadius())) {
                    return distanceEstimate(sub, point, k, comp);
                }
            }
            // okay, this directory node contains the point, but none of the leaves does.
            // We just return the size of this node, because all it's leaf nodes should
            // contain more than enough candidate in proximity of 'point'.
            return node->getRadius() * std::sqrt(point.size()); // TODO backport???  sqrt(3) ?!?!
        }
      }

    private:
    double rangeSearchKNN(
        const QNode<Key, T>* node,
        const Key& center,
        std::vector<QEntryDist<Key, T>>& candidates,
        size_t k,
        double maxRange) {
        if (node->isLeaf()) {
            auto& points = node->getEntries();
            for (size_t i = 0; i < points.size(); i++) {
                QEntry<Key, T>& p = points[i];
                double dist = distance(center, p.point());
                if (dist < maxRange) {
                    candidates.add(new QEntryDist<Key, T>(p, dist));
                }
            }
            maxRange = adjustRegionKNN(candidates, k, maxRange);
        } else {
            std::vector<QNode<Key, T>>& nodes = node->getChildNodes();
            for (size_t i = 0; i < nodes.size(); i++) {
                QNode<Key, T>* sub = nodes[i];
                if (sub != nullptr &&
                    distToRectNode(center, sub->getCenter(), sub->getRadius()) < maxRange) {
                    maxRange = rangeSearchKNN(sub, center, candidates, k, maxRange);
                    // we set maxRange simply to the latest returned value.
                }
            }
        }
        return maxRange;
    }

  private:
    double adjustRegionKNN(std::vector<QEntryDist<Key, T>>& candidates, size_t k, double maxRange) {
        if (candidates.size() < k) {
            // wait for more candidates
            return maxRange;
        }

        // use stored distances instead of recalculating them
        auto comp = [](const QEntryDist<Key, T>& c1, const QEntryDist<Key, T>& c2) {
            return c1.dist() < c2.dist();
        };
        std::sort(candidates.begin(), candidates.end(), comp); // TODO why are we sorting the whole list???
        //candidates.sort(QEntryDist.COMP);

        while (candidates.size() > k) {
            candidates.erase(candidates.end() - 1);
        }

        double range = candidates.back().dist();
        return range;
    }

//  private:
//    class QQueryIteratorKNN { //implements QueryIteratorKNN<PointEntryDist<T>> {
//        using IterT = decltype(std::vector<QEntryDist<Key, T>>().begin());
//        IterT it;
//
//      public:
//        QQueryIteratorKNN(const Key& center, int k) {
//            it = knnQuery(center, k).begin();
//        }
//
//        public: bool hasNext() {
//            return it.hasNext();
//        }
//
//        public: const QEntryDist<Key, T>& next() {
//            return it.next();
//        }
//    };

    /**
	 * Returns a printable list of the tree.
	 * @return the tree as String
     */
//    public: String
//    toStringTree() {
//        StringBuilder sb = new StringBuilder();
//        if (root_ == nullptr) {
//            sb.append("empty tree");
//        } else {
//            toStringTree(sb, root, 0, 0);
//        }
//        return sb.toString();
//    }
//
//     private: void toStringTree(
//        StringBuilder sb, QNode<Key, T>* node, int depth, int posInParent) {
//        Iterator < ? > it = node->getChildIterator();
//        String prefix = "";
//        for (size_t i = 0; i < depth; i++) {
//            prefix += ".";
//        }
//        sb.append(prefix + posInParent + " d=" + depth);
//        sb.append(" " + Arrays.toString(node->getCenter()));
//        sb.append("/" + node->getRadius() + NL);
//        prefix += " ";
//        int pos = 0;
//        while (it.hasNext()) {
//            Object o = it.next();
//            if (o instanceof QNode) {
//                QNode<Key, T>* sub = (QNode<Key, T>)o;
//                toStringTree(sb, sub, depth + 1, pos);
//            } else if (o instanceof QEntry) {
//                QEntry<Key, T>& e = (QEntry<Key, T>)o;
//                sb.append(prefix + Arrays.toString(e.point()));
//                sb.append(" v=" + e.value() + NL);
//            }
//            pos++;
//        }
//    }
//
//    public: String toString() {
//        return "QuadTreeKD0;maxNodeSize=" + maxNodeSize + ";maxDepth=" + MAX_DEPTH +
//            ";DEBUG=" + DEBUG + ";center/radius=" +
//            (root_ == nullptr ? "nullptr" : (Arrays.toString(root_->getCenter()) + "/" + root_->getRadius()));
//    }

    public: QStats getStats() {
        QStats s{};
        if (root_ != nullptr) {
            root_->checkNode(s, nullptr, 0);
        }
        return s;
    }

    public: int
    getDims() {
        return dims;
    }

    public:
      auto begin() {
        if (root_ == nullptr) {
            return begin_query({Key{}, Key{}});
        }
        // return query(root_->);
        // TODO
        assert(false);
    }

    bool empty() {
        return size_ == 0;
    }

    template <typename CALLBACK, typename FILTER = FilterNoOp>
    void for_each(QueryBox query_box, CALLBACK&& callback, FILTER&& filter = FILTER()) const {
        auto it = begin_query(query_box, std::forward<FILTER>(filter));
        // TODO move FILTER into begin_query
        auto end = this->end();
        while (it != end) {
            Key k = it._entry()->point();
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
        return QIterator<Key, T, FILTER>(
            root_, query_box.min(), query_box.max(), std::forward<FILTER>(filter));
    }

    auto end() const {
        return QIteratorFind<Key, T>(nullptr, nullptr);
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
        return QIteratorKnn<Key, T>(std::move(result));
    }


//    public: QQueryIteratorKNN queryKNN(const Key& center, int k) {
//        return new QQueryIteratorKNN(center, k);
//    }

    public: int getNodeCount() {
        return getStats().getNodeCount();
    }

    public: int getDepth() {
        return getStats().getMaxDepth();
    }
};

}
#endif  // TINSPIN_QUADTREE_POINT_H
