// SPDX-FileCopyrightText: 2023 Tilmann Zäschke <zoodb@gmx.de>
// SPDX-License-Identifier: MIT

#ifndef TINSPIN_QUADTREE_POINT_H
#define TINSPIN_QUADTREE_POINT_H

#include "include/phtree/common/common.h"
#include "include/phtree/converter.h"
#include "include/phtree/filter.h"
#include "src/util/ph-util.h"
#include <cmath>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_set>

namespace tinspin {

using namespace improbable::phtree;

/*
 * A simple MX-quadtree implementation with configurable maximum depth, maximum nodes size, and
 * (if desired) automatic guessing of root rectangle.
 *
 * This version of the quadtree stores for each node only the center point and the
 * distance (radius) to the edges.
 * This reduces space requirements but increases problems with numerical precision.
 * Overall it is more space efficient and slightly faster.
 */

namespace {

template <typename KEy, typename T>
class QREntry;

static const double EPS_MUL = 1.000000001;

template <typename Key>
static bool isPointEnclosed(const Key& point, const Key& min, const Key& max) {
    for (size_t d = 0; d < min.size(); ++d) {
        if (point[d] < min[d] || point[d] > max[d]) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isPointEnclosed(const Key& point, const Key& center, double radius) {
    for (size_t d = 0; d < center.size(); ++d) {
        if (point[d] < center[d] - radius || point[d] > center[d] + radius) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isPointEqual(const Key& p1, const Key& p2) {
    for (size_t d = 0; d < p1.size(); ++d) {
        if (p1[d] != p2[d]) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isRectEqual(const Key& p1L, const Key& p1U, const Key& p2L, const Key& p2U) {
    return isPointEqual(p1L, p2L) && isPointEqual(p1U, p2U);
}

template <typename Key, typename T>
static bool isRectEqual(QREntry<Key, T> e, const Key& keyL, const Key& keyU) {
    return isRectEqual(e.lower(), e.upper(), keyL, keyU);
}

template <typename Key>
static bool overlap(const Key& min, const Key& max, const Key& min2, const Key& max2) {
    for (size_t d = 0; d < min.size(); ++d) {
        if (max[d] < min2[d] || min[d] > max2[d]) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool overlap(const Key& min, const Key& max, const Key& center, double radius) {
    for (size_t d = 0; d < min.size(); ++d) {
        if (max[d] < center[d] - radius || min[d] > center[d] + radius) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool overlap(const Key& center, double radius, const Key& center2, double radius2) {
    for (size_t d = 0; d < center.size(); ++d) {
        if (center[d] + radius < center2[d] - radius2 ||
            center[d] - radius > center2[d] + radius2) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isRectEnclosed(
    const Key& minEnclosed, const Key& maxEnclosed, const Key& minOuter, const Key& maxOuter) {
    for (size_t d = 0; d < minOuter.size(); ++d) {
        if (maxOuter[d] < maxEnclosed[d] || minOuter[d] > minEnclosed[d]) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isRectEnclosed(
    const Key& minEnclosed, const Key& maxEnclosed, const Key& centerOuter, double radiusOuter) {
    for (size_t d = 0; d < centerOuter.size(); ++d) {
        double radOuter = radiusOuter;
        if ((centerOuter[d] + radOuter) < maxEnclosed[d] ||
            (centerOuter[d] - radOuter) > minEnclosed[d]) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isRectEnclosed(
    const Key& centerEnclosed, double radiusEnclosed, const Key& centerOuter, double radiusOuter) {
    for (size_t d = 0; d < centerOuter.size(); ++d) {
        double radOuter = radiusOuter;
        double radEncl = radiusEnclosed;
        if ((centerOuter[d] + radOuter) < (centerEnclosed[d] + radEncl) ||
            (centerOuter[d] - radOuter) > (centerEnclosed[d] - radEncl)) {
            return false;
        }
    }
    return true;
}

template <typename Key>
static bool isRectEnclosed(
    const Key& centerEnclosed, double radiusEnclosed, const Key& minOuter, const Key& maxOuter) {
    for (size_t d = 0; d < centerEnclosed.size(); ++d) {
        double radEncl = radiusEnclosed;
        if (maxOuter[d] < (centerEnclosed[d] + radEncl) ||
            minOuter[d] > (centerEnclosed[d] - radEncl)) {
            return false;
        }
    }
    return true;
}

// /**
// * Calculates distance to center point of rectangle.
// * @param p point
// * @param rMin rectangle min
// * @param rMax rectangle max
// * @return distance to center point
// */
// template <typename Key>
// static double distToRectCenter(const Key& p, const Key& rMin, const Key& rMax) {
//    double dist = 0;
//    for (size_t i = 0; i < p.size(); ++i) {
//        double d = (rMin[i] + rMax[i]) / 2. - p[i];
//        dist += d * d;
//    }
//    return std::sqrt(dist);
//}
//
// /**
// * Calculates distance to center point of rectangle.
// * @param p point
// * @param e rectangle
// * @return distance to center point
// */
// template <typename Key, typename T>
// static double distToRectCenter(const Key& p, QREntry<Key, T> e) {
//    return distToRectCenter(p, e.lower(), e.upper());
//}
//
// /**
// * Calculates distance to the edge of rectangle.
// * @param p point
// * @param rMin rectangle min
// * @param rMax rectangle max
// * @return distance to edge
// */
// template <typename Key>
// double distToRectEdge(const Key& center, const Key& rLower, const Key& rUpper) {
//    double dist = 0;
//    for (size_t i = 0; i < center.size(); ++i) {
//        double d = 0;
//        if (center[i] > rUpper[i]) {
//            d = center[i] - rUpper[i];
//        } else if (center[i] < rLower[i]) {
//            d = rLower[i] - center[i];
//        }
//        dist += d * d;
//    }
//    return std::sqrt(dist);
//}
//
///**
// * Calculates distance to edge of rectangle.
// * @param p point
// * @param e rectangle
// * @return distance to edge point
// */
// template <typename Key, typename T>
// static double distToRectEdge(const Key& p, QREntry<Key, T> e) {
//    return distToRectEdge(p, e.lower(), e.upper());
//}

/**
 * Calculates distance to the edge of a node.
 * @param point the point
 * @param nodeCenter the center of the node
 * @param nodeRadius radius of the node
 * @return distance to edge of the node or 0 if the point is inside the node
 */
template <typename Key, typename DISTANCE>
static double distToRectNode(
    const Key& point, const Key& nodeCenter, double nodeRadius, DISTANCE& dist_fn) {
    Key closest{point};
    for (size_t i = 0; i < point.size(); ++i) {
        if (point[i] > nodeCenter[i] + nodeRadius) {
            closest[i] = nodeCenter[i] + nodeRadius;
        } else if (point[i] < nodeCenter[i] - nodeRadius) {
            // TODO precision??
            closest[i] = nodeCenter[i] - nodeRadius;
        }
    }
    return dist_fn(point, closest);
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
    Key point_;
    T value_;

  public:
    QEntry(const Key& key, const T& value) : point_{key}, value_{value} {}

    const Key& key() const {
        return point_;
    }

    const T& value() const {
        return value_;
    }

    void setKey(const Key& newPoint) {
        point_ = newPoint;
    }
};

template <typename Key, typename T>
class QEntryDist {
  public:
    QEntryDist(QEntry<Key, T>* e, double dist) : entry_{e}, distance_{dist} {}

    const Key& key() const {
        return entry_->key();
    }

    const T& value() const {
        return entry_->value();
    }

    double dist() const {
        return distance_;
    }

    QEntry<Key, T>* _entry() {
        return entry_;
    }

  private:
    QEntry<Key, T>* entry_;
    double distance_;
};

/**
 * Node class for the quadtree.
 */
template <typename Key, typename T>
class QNode {
    Key center_;
    double radius_;
    // null indicates that we have sub-nopde i.o. values
    std::vector<QEntry<Key, T>> values_;
    std::vector<QNode<Key, T>*> subs_;
    bool is_leaf_;

  public:
    QNode(const Key& center, double radius)
    : center_{center}, radius_{radius}, values_{}, subs_{}, is_leaf_{true} {
        values_.reserve(2);
    }

    QNode(const Key& center, double radius, QNode<Key, T>* subNode)
    : center_{center}, radius_{radius}, values_{}, subs_{}, is_leaf_{false} {
        subs_.emplace_back(subNode);
    }

    ~QNode() {
        for (auto* sub : subs_) {
            delete sub;
        }
    }

    QNode<Key, T>* tryPut(const Key& key, const T& value, size_t maxNodeSize, bool enforceLeaf) {
        assert(isPointEnclosed(key, center_, radius_));

        // traverse subs?
        if (!isLeaf()) {
            return getOrCreateSub(key);
        }

        // add if:
        // a) we have space
        // b) we have maxDepth
        // c) elements are equal (work only for n=1, avoids splitting
        //    in cases where splitting won't help. For n>1 the
        //    local limit is (temporarily) violated. // TODO we should check all points!?
        if (values_.size() < maxNodeSize || enforceLeaf || key == values_[0].key()) {
            values_.emplace_back(key, value);
            return nullptr;
        }

        assert(subs_.empty());
        for (size_t i = 0; i < values_.size(); ++i) {
            QEntry<Key, T>& e2 = values_[i];
            QNode<Key, T>* sub = getOrCreateSub(e2.key());
            while (sub != nullptr) {
                // This may recurse if all entries fall
                // into the same subnode
                // TODO std::move key/value?
                sub = sub->tryPut(e2.key(), e2.value(), maxNodeSize, false);
            }
        }
        values_.clear();
        values_.shrink_to_fit();
        is_leaf_ = false;
        return getOrCreateSub(key);
    }

  private:
    QNode<Key, T>* getOrCreateSub(const Key& key) {
        QNode<Key, T>* n = findSubNode(key);
        if (n == nullptr) {
            n = createSubForEntry(key);
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
        for (size_t i = 0; i < subs_.size(); ++i) {
            QNode<Key, T>* n = subs_[i];
            if (isPointEnclosed(p, n->center_, n->radius_)) {
                return n;
            }
        }
        return nullptr;
    }

  public:
    size_t remove(QNode<Key, T>* parent, const Key& key, size_t maxNodeSize) {
        if (!is_leaf_) {
            QNode<Key, T>* sub = findSubNode(key);
            if (sub != nullptr) {
                return sub->remove(this, key, maxNodeSize);
            }
            return 0;
        }

        size_t n = 0;
        for (size_t i = 0; i < values_.size(); ++i) {
            QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.key(), key)) {
                values_.erase(values_.begin() + i);
                ++n;
                --i;
            }
        }
        // TODO provide threshold for re-insert
        // i.e. do not always merge.
        if (n > 0 && parent != nullptr) {
            parent->checkAndMergeLeafNodes(maxNodeSize);
        }
        return n;
    }

    size_t remove(QNode<Key, T>* parent, const Key& key, size_t maxNodeSize, const T& value) {
        if (!is_leaf_) {
            QNode<Key, T>* sub = findSubNode(key);
            if (sub != nullptr) {
                return sub->remove(this, key, maxNodeSize, value);
            }
            return 0;
        }

        size_t n = 0;
        for (size_t i = 0; i < values_.size(); ++i) {
            QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.key(), key) && e.value() == value) {
                values_.erase(values_.begin() + i);
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

    /**
     * @return 1 == found
     */
    size_t update(
        QNode<Key, T>* parent,
        const Key& keyOld,
        const Key& keyNew,
        const T& value,
        size_t maxNodeSize,
        bool& requiresReinsert,
        size_t currentDepth,
        size_t maxDepth) {
        if (!is_leaf_) {
            QNode<Key, T>* sub = findSubNode(keyOld);
            if (sub == nullptr) {
                return 0;
            }
            size_t ret = sub->update(
                this,
                keyOld,
                keyNew,
                value,
                maxNodeSize,
                requiresReinsert,
                currentDepth + 1,
                maxDepth);
            if (ret != 0 && requiresReinsert && isPointEnclosed(keyNew, center_, radius_)) {
                requiresReinsert = false;
                auto* r = this;  // TODO backport, r is always a QNode
                while (r != nullptr) {
                    r = r->tryPut(keyNew, value, maxNodeSize, currentDepth++ > maxDepth);
                }
            }
            return ret;  // TODO backport -> return 0 if entry was not found!
        }

        for (size_t i = 0; i < values_.size(); ++i) {
            QEntry<Key, T>& e = values_[i];
            if (isPointEqual(e.key(), keyOld) && e.value() == value) {
                if (isPointEnclosed(keyNew, center_, radius_)) {
                    // reinsert locally;
                    // TODO avoid previous deletion / or at least use std::move!!!!!!!
                    //    --> backport
                    // values_.emplace_back(e);
                    e.setKey(keyNew);
                    requiresReinsert = false;
                } else {
                    values_.erase(values_.begin() + i);
                    requiresReinsert = true;
                    // TODO provide threshold for re-insert
                    // i.e. do not always merge.
                    if (parent != nullptr) {
                        parent->checkAndMergeLeafNodes(maxNodeSize);
                    }
                }
                return 1;
            }
        }
        requiresReinsert = false;
        return 0;
    }

  private:
    void checkAndMergeLeafNodes(size_t maxNodeSize) {
        // check
        size_t nTotal = 0;
        for (size_t i = 0; i < subs_.size(); ++i) {
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
        // values_ = new ArrayList<>(nTotal);
        for (size_t i = 0; i < subs_.size(); ++i) {
            values_.insert(values_.begin(), subs_[i]->values_.begin(), subs_[i]->values_.end());
        }
        // subs = nullptr;
        for (auto sub : subs_) {
            delete sub;
        }
        subs_.clear();
        subs_.shrink_to_fit();
        is_leaf_ = true;
    }

  public:
    const Key& getCenter() const noexcept {
        return center_;
    }

    double getRadius() const noexcept {
        return radius_;
    }

    const QNode<Key, T>* getExactLeaf(const Key& key) const {
        if (is_leaf_) {
            return this;
        }
        QNode<Key, T>* sub = findSubNode(key);
        if (sub != nullptr) {
            return sub->getExactLeaf(key);
        }
        return nullptr;
    }

    auto& entries() {
        return values_;
    }

    auto& entries() const {
        return values_;
    }

    void stats(QStats s, QNode<Key, T>* parent, size_t depth) {
        if (depth > s.maxDepth) {
            s.maxDepth = depth;
        }
        s.nNodes++;

        if (parent != nullptr) {
            if (!isRectEnclosed(center_, radius_, parent->center_, parent->radius * EPS_MUL)) {
                for (size_t d = 0; d < center_.size(); ++d) {
                    //					if ((centerOuter[d]+radiusOuter) /
                    //(centerEnclosed[d]+radiusEnclosed) < 0.9999999 ||
                    //(centerOuter[d]-radiusOuter) / (centerEnclosed[d]-radiusEnclosed) > 1.0000001)
                    //{ 						return false;
                    //					}
                    std::cout << "Outer: " << parent->radius_ << " " << parent->center_
                              << std::endl;
                    std::cout << "Child: " << radius_ << " " << center_ << std::endl;
                    std::cout << (parent->center_[d] + parent->radius_) << " vs "
                              << (center_[d] + radius_) << std::endl;
                    std::cout << "r="
                              << (parent->center_[d] + parent->radius_) / (center_[d] + radius_)
                              << std::endl;
                    std::cout << (parent->center_[d] - parent->radius_) << " vs "
                              << (center_[d] - radius_) << std::endl;
                    std::cout << "r="
                              << (parent->center_[d] - parent->radius_) / (center_[d] - radius_)
                              << std::endl;
                }
                assert(false);
            }
        }
        if (values_ != nullptr) {
            for (size_t i = 0; i < values_.size(); ++i) {
                QEntry<Key, T>& e = values_[i];
                if (!isPointEnclosed(e.key(), center_, radius_ * EPS_MUL)) {
                    std::cout << "Node: " << radius_ << " " << center_ << std::endl;
                    std::cout << "Child: " << e.key() << std::endl;
                    for (size_t d = 0; d < center_.size(); ++d) {
                        //						if ((centerOuter[d]+radiusOuter) /
                        //(centerEnclosed[d]+radiusEnclosed) < 0.9999999 ||
                        //								(centerOuter[d]-radiusOuter) /
                        //(centerEnclosed[d]-radiusEnclosed) > 1.0000001) {
                        // return false;
                        //						}
                        std::cout << "min/max for " << d << std::endl;
                        std::cout << "min: " << (center_[d] - radius_) << " vs " << (e.key()[d])
                                  << std::endl;
                        std::cout << "r=" << (center_[d] - radius_) / (e.key()[d]) << std::endl;
                        std::cout << "max: " << (center_[d] << radius_) << " vs " << (e.key()[d])
                                  << std::endl;
                        std::cout << "r=" << (center_[d] << radius_) / (e.key()[d]) << std::endl;
                    }
                    assert(false);
                }
            }
            if (subs_ != nullptr) {
                assert(false);
            }
        } else {
            for (size_t i = 0; i < subs_.size(); ++i) {
                QNode<Key, T>& n = subs_[i];
                n.stats(s, this, depth + 1);
            }
        }
    }

    bool isLeaf() const noexcept {
        return is_leaf_;
    }

    auto& getChildNodes() {
        assert(!isLeaf());
        return subs_;
    }

    auto& getChildNodes() const {
        assert(!isLeaf());
        return subs_;
    }
};

template <typename Key, typename T>
class QIteratorBase {
  public:
    explicit QIteratorBase() noexcept : entry_{nullptr} {}

    inline auto& operator*() const noexcept {
        assert(!IsEnd() && "This iterator is invalid.");
        return entry_->value();
    }

    inline auto* operator->() const noexcept {
        assert(!IsEnd() && "This iterator is invalid.");
        return &entry_->value();
    }

    inline friend bool operator==(
        const QIteratorBase<Key, T>& left, const QIteratorBase<Key, T>& right) noexcept {
        return left.entry_ == right.entry_;
    }

    inline friend bool operator!=(
        const QIteratorBase<Key, T>& left, const QIteratorBase<Key, T>& right) noexcept {
        return left.entry_ != right.entry_;
    }

    auto _entry() const noexcept {
        return entry_;
    }

  protected:
    bool IsEnd() const noexcept {
        return entry_ == nullptr;
    }

    void SetFinished() noexcept {
        entry_ = nullptr;
    }

    void SetCurrentResult(QEntry<Key, T>* entry) noexcept {
        entry_ = entry;
    }

  protected:
    QEntry<Key, T>* entry_ = nullptr;
};

template <typename Key, typename T, typename FILTER>
class QIterator : public QIteratorBase<Key, T> {
    using IterNodeT = decltype(std::vector<QNode<Key, T>*>().cbegin());
    using IterEntryT = decltype(std::vector<QEntry<Key, T>>().cbegin());
    using EntryInnerT = std::pair<QNode<Key, T>*, IterNodeT>;

    std::stack<EntryInnerT> stack_;  // TODO backport, Why is this a Deque????????
    IterEntryT iter_leaf_;
    Key min;
    Key max;
    FILTER filter_;

  public:
    template <typename F = FilterNoOp>
    QIterator(QNode<Key, T>* root, const Key& min, const Key& max, F&& filter = F())
    : QIteratorBase<Key, T>()
    , stack_{}
    , min(min)
    , max(max)
    , filter_(std::forward<F>(filter)) {
        if (root != nullptr) {
            if (root->isLeaf()) {
                iter_leaf_ = root->entries().begin();
                stack_.emplace(std::make_pair(root, IterNodeT{}));
            } else {
                stack_.emplace(std::make_pair(root, root->getChildNodes().begin()));
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

  private:
    void findNext() {
        while (!stack_.empty()) {
            // Are we currently iterating a leaf?
            if (stack_.top().first->isLeaf()) {
                if (findNextInNode()) {
                    return;
                }
                stack_.pop();
                if (stack_.empty()) {
                    this->SetFinished();
                    return;
                }
                ++stack_.top().second;
            }

            // traverse inner nodes
            auto& ee = stack_.top();
            auto& it = ee.second;
            while (it != ee.first->getChildNodes().end()) {
                auto* node = *it;
                if (overlap(min, max, node->getCenter(), node->getRadius())) {
                    if (node->isLeaf()) {
                        iter_leaf_ = node->entries().begin();
                        stack_.push(std::make_pair(node, IterNodeT{}));
                        if (findNextInNode()) {
                            return;
                        }
                        stack_.pop();
                    } else {
                        IterNodeT it2 = node->getChildNodes().begin();
                        stack_.push(std::make_pair(node, it2));
                        break;
                    }
                }
                ++it;
            }
            if (it == ee.first->getChildNodes().end()) {
                stack_.pop();
                if (!stack_.empty()) {
                    ++stack_.top().second;
                }
            }
        }
        this->SetFinished();
    }

    bool findNextInNode() {
        auto* node = stack_.top().first;
        assert(node->isLeaf());
        while (iter_leaf_ != node->entries().end()) {
            const QEntry<Key, T>& e = *iter_leaf_;
            ++iter_leaf_;
            if (isPointEnclosed(e.key(), min, max) && filter_.IsEntryValid(e.key(), e.value())) {
                this->SetCurrentResult(const_cast<QEntry<Key, T>*>(&e));
                return true;
            }
        }
        return false;
    }
};

template <typename Key, typename T, typename FILTER>
class QIteratorFind : public QIteratorBase<Key, T> {
    using IterEntryT = decltype(std::vector<QEntry<Key, T>>().begin());

    IterEntryT iter_;
    QNode<Key, T>* leaf_;
    FILTER filter_;

  public:
    template <typename F>
    QIteratorFind(QNode<Key, T>* leaf, F&& filter = F())
    : QIteratorBase<Key, T>(), leaf_{leaf}, filter_(std::forward<F>(filter)) {
        if (leaf != nullptr) {
            iter_ = leaf->entries().begin();
            findNext();
        } else {
            this->SetFinished();
        }
    }

    QIteratorFind& operator++() noexcept {
        assert(!this->IsEnd());
        findNext();
        return *this;
    }

    QIteratorFind operator++(int) noexcept {
        assert(!this->IsEnd());
        QIterator iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void findNext() {
        while (iter_ != leaf_->entries().end()) {
            auto& entry = *iter_;
            ++iter_;
            if (filter_(entry)) {
                this->SetCurrentResult(&entry);
                return;
            }
        }
        this->SetFinished();
    }
};

template <typename Key, typename T>
class QIteratorEnd : public QIteratorBase<Key, T> {
  public:
    template <typename F = FilterNoOp>
    QIteratorEnd() : QIteratorBase<Key, T>() {
        this->SetFinished();
    }

    QIteratorEnd& operator++() noexcept {
        assert(false);
        return *this;
    }

    QIteratorEnd operator++(int) noexcept {
        assert(false);
        QIterator iterator(*this);
        ++(*this);
        return iterator;
    }
};

template <typename Key, typename T, typename DISTANCE, typename FILTER>
class QIteratorKnnHS : public QIteratorBase<Key, T> {
    static constexpr dimension_t DIM = 3;
    struct EntryDistT {
        double first;                  // distance
        const QEntry<Key, T>* second;  // entry
        const QNode<Key, T>* node;     // node

        EntryDistT(double dist, const QEntry<Key, T>* e) : first{dist}, second{e}, node{nullptr} {}
        EntryDistT(double dist, const QNode<Key, T>* node)
        : first{dist}, second{nullptr}, node{node} {}

        [[nodiscard]] bool is_node() const {
            return node != nullptr;
        }
    };

    struct CompareEntryDistByDistance {
        bool operator()(const EntryDistT& left, const EntryDistT& right) const {
            return left.first > right.first;
        };
    };

  public:
    template <typename DIST, typename F>
    explicit QIteratorKnnHS(
        const QNode<Key, T>* root,
        size_t min_results,
        const Key& center,
        DIST&& dist_fn,
        F&& filter_fn)
    : QIteratorBase<Key, T>()
    , center_{center}
    , current_distance_{std::numeric_limits<double>::max()}
    , num_found_results_(0)
    , num_requested_results_(min_results)
    , filter_(std::forward<F>(filter_fn))
    , distance_(std::forward<DIST>(dist_fn)) {
        if (min_results <= 0 || root == nullptr) {
            this->SetFinished();
            return;
        }

        double dist = distToRectNode(center, root->getCenter(), root->getRadius(), distance_);
        queue_.emplace(dist, root);
        FindNextElement();
    }

    [[nodiscard]] double distance() const {
        return current_distance_;
    }

    const Key& first() const noexcept {
        return this->entry_->key();
    }

    QIteratorKnnHS& operator++() noexcept {
        FindNextElement();
        return *this;
    }

    QIteratorKnnHS operator++(int) noexcept {
        QIteratorKnnHS iterator(*this);
        ++(*this);
        return iterator;
    }

  private:
    void FindNextElement() {
        while (num_found_results_ < num_requested_results_ && !queue_.empty()) {
            auto& candidate = queue_.top();
            if (!candidate.is_node()) {
                // data entry
                ++num_found_results_;
                this->SetCurrentResult(const_cast<QEntry<Key, T>*>(candidate.second));
                current_distance_ = candidate.first;
                // We need to pop() AFTER we processed the value, otherwise the reference is
                // overwritten.
                queue_.pop();
                return;
            } else {
                // inner node
                auto* node = candidate.node;
                queue_.pop();
                if (node->isLeaf()) {
                    for (auto& entry : node->entries()) {
                        if (filter_.IsEntryValid(entry.key(), entry.value())) {
                            double d = distance_(center_, entry.key());
                            queue_.emplace(d, &entry);
                        }
                    }
                } else {
                    for (auto* subnode : node->getChildNodes()) {
                        double dist = distToRectNode(
                            center_, subnode->getCenter(), subnode->getRadius(), distance_);
                        queue_.emplace(dist, subnode);
                    }
                }
            }
        }
        this->SetFinished();
        current_distance_ = std::numeric_limits<double>::max();
    }

  private:
    const Key center_;
    double current_distance_;
    std::priority_queue<EntryDistT, std::vector<EntryDistT>, CompareEntryDistByDistance> queue_;
    size_t num_found_results_;
    size_t num_requested_results_;
    FILTER filter_;
    DISTANCE distance_;
};

template <typename Key, typename T, typename CALLBACK, typename FILTER>
class ForEach {
  public:
    template <typename CB, typename F>
    ForEach(const Key& range_min, const Key& range_max, CB&& callback, F&& filter)
    : range_min_{range_min}
    , range_max_{range_max}
    , callback_{std::forward<CB>(callback)}
    , filter_(std::forward<F>(filter)) {}

    void Traverse(const QNode<Key, T>* parent_node) {
        if (parent_node != nullptr) {
            if (parent_node->isLeaf()) {
                TraverseLeaf(parent_node, false);

            } else {
                TraverseInner(parent_node, false);
            }
        }
    }

  private:
    void TraverseInner(const QNode<Key, T>* parent_node, const bool is_enclosed) {
        for (const auto node : parent_node->getChildNodes()) {
            if (is_enclosed ||
                overlap(range_min_, range_max_, node->getCenter(), node->getRadius())) {
                bool is_node_enclosed = is_enclosed ||
                    isRectEnclosed(node->getCenter(), node->getRadius(), range_min_, range_max_);
                if (!node->isLeaf()) {
                    TraverseInner(node, is_node_enclosed);
                } else {
                    TraverseLeaf(node, is_node_enclosed);
                }
            }
        }
    }

    void TraverseLeaf(const QNode<Key, T>* node, const bool is_enclosed) {
        for (auto& entry : node->entries()) {
            const Key& key = entry.key();
            const T& value = entry.value();
            // TODO skip if node is fully enclosed!
            if ((is_enclosed || isPointEnclosed(key, range_min_, range_max_)) &&
                filter_.IsEntryValid(key, value)) {
                callback_(key, value);
            }
        }
    }

    const Key range_min_;
    const Key range_max_;
    CALLBACK callback_;
    FILTER filter_;
};

}  // namespace

template <typename Key, typename T>
class QuadTree {
    static const size_t MAX_DEPTH = 50;
    using QueryBox = tinspin::Box<Key>;

    static const size_t DEFAULT_MAX_NODE_SIZE = 10;

    const size_t dims;
    const size_t maxNodeSize;
    QNode<Key, T>* root_ = nullptr;
    size_t size_ = 0;

  public:
    using KeyInternal = Key;

    QuadTree(size_t dims = 3, size_t maxNodeSize = DEFAULT_MAX_NODE_SIZE)
    : dims{dims}, maxNodeSize{maxNodeSize} {}

    ~QuadTree() {
        delete root_;
    }

    /**
     * Insert a key-value pair.
     * @param key the key
     * @param value the value
     */
    template <typename T2>
    void emplace(const Key& key, T2&& value) {
        size_++;
        if (root_ == nullptr) {
            initializeRoot(key);
        }
        ensureCoverage(key);
        auto* r = root_;
        size_t depth = 0;
        while (r != nullptr) {  // TODO backport, r is always a QNode
            r = r->tryPut(key, std::forward<T2>(value), maxNodeSize, depth++ > MAX_DEPTH);
        }
    }

    void insert(const Key& key, const T& value) {
        size_++;
        if (root_ == nullptr) {
            initializeRoot(key);
        }
        ensureCoverage(key);
        auto* r = root_;
        size_t depth = 0;
        while (r != nullptr) {  // TODO backport, r is always a QNode
            r = r->tryPut(key, value, maxNodeSize, depth++ > MAX_DEPTH);
        }
    }

  private:
    void initializeRoot(const Key& key) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (size_t d = 0; d < dims; ++d) {
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
        for (size_t d = 0; d < dims; ++d) {
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
        size_t n = 0;
        for_each({key, key}, [&n](const Key&, const T&) { ++n; });
        return n;
    }

    /**
     * Get the value associates with the key.
     * @param key the key to look up
     * @return the value for the key or 'nullptr' if the key was not found
     */
    auto find(const Key& key) const {
        auto filter = [&key](const QEntry<Key, T>& e) { return key == e.key(); };
        if (root_ == nullptr) {
            return QIteratorFind<Key, T, decltype(filter)>(nullptr, filter);
        }
        QNode<Key, T>* leaf = const_cast<QNode<Key, T>*>(root_->getExactLeaf(key));
        return QIteratorFind<Key, T, decltype(filter)>(leaf, filter);
    }

    auto find(const Key& key, const T& value) const {
        auto filter = [&key, &value](const QEntry<Key, T>& e) {
            return key == e.key() && value == e.value();
        };
        if (root_ == nullptr) {
            return QIteratorFind<Key, T, decltype(filter)>(nullptr, filter);
        }
        QNode<Key, T>* leaf = const_cast<QNode<Key, T>*>(root_->getExactLeaf(key));
        return QIteratorFind<Key, T, decltype(filter)>(leaf, filter);
    }

    /**
     * Remove a key.
     * @param key key to remove
     * @return the value associated with the key or 'nullptr' if the key was not found
     */
    size_t erase(const Key& key) {
        if (root_ == nullptr) {
            return 0;
        }
        size_t n = root_->remove(nullptr, key, maxNodeSize);
        if (n == 0) {
            return 0;
        }
        size_ -= n;
        return n;
    }

    size_t erase(const Key& key, const T& value) {
        if (root_ == nullptr) {
            return 0;
        }
        size_t n = root_->remove(nullptr, key, maxNodeSize, value);
        if (n == 0) {
            return 0;
        }
        size_ -= n;
        return n;
    }

    /**
     * Reinsert the key.
     * @param oldKey old key
     * @param newKey new key
     * @return the value associated with the key or 'nullptr' if the key was not found.
     */
    size_t relocate(const Key& oldKey, const Key& newKey, const T& value) {
        if (root_ == nullptr) {
            return 0;
        }
        bool requiresReinsert = false;
        size_t result = root_->update(
            nullptr, oldKey, newKey, value, maxNodeSize, requiresReinsert, 0, MAX_DEPTH);
        if (result == 0) {
            // not found
            return 0;
        }
        if (requiresReinsert) {
            // does not fit in root node...
            ensureCoverage(newKey);
            auto* r = root_;
            size_t depth = 0;
            while (r != nullptr) {  // TODO backport, r is always a QNode
                r = r->tryPut(newKey, value, maxNodeSize, depth++ > MAX_DEPTH);
            }
        }
        return 1;
    }

    /**
     * Ensure that the tree covers the entry.
     * @param e Entry to cover.
     */
  private:
    // TODO backport : pass in Key only.
    void ensureCoverage(const Key& p) {
        while (!isPointEnclosed(p, root_->getCenter(), root_->getRadius())) {
            const Key& center = root_->getCenter();
            double radius = root_->getRadius();
            Key center2{};
            double radius2 = radius * 2;
            // TODO use DIM?
            for (size_t d = 0; d < center.size(); ++d) {
                if (p[d] < center[d] - radius) {
                    center2[d] = center[d] - radius;
                    // root will end up in upper quadrant in this
                    // dimension
                } else {
                    // extend upwards, even if extension unnecessary for this dimension.
                    center2[d] = center[d] + radius;
                }
            }
            assert(isRectEnclosed(center, radius, center2, radius2));
            root_ = new QNode<Key, T>(center2, radius2, root_);
        }
    }

    /**
     * Get the number of key-value pairs in the tree.
     * @return the size
     */
  public:
    size_t size() {
        return size_;
    }

    /**
     * Removes all elements from the tree.
     */
    void clear() {
        size_ = 0;
        delete root_;
        root_ = nullptr;
    }

    QStats stats() {
        QStats s{};
        if (root_ != nullptr) {
            root_->stats(s, nullptr, 0);
        }
        return s;
    }

    size_t dimensions() {
        return dims;
    }

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
        ForEach<Key, T, CALLBACK, FILTER>(
            query_box.min(),
            query_box.max(),
            std::forward<CALLBACK>(callback),
            std::forward<FILTER>(filter))
            .Traverse(root_);
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
        return QIteratorEnd<Key, T>();
    }

    template <typename DISTANCE, typename FILTER = FilterNoOp>
    auto begin_knn_query(
        size_t k,
        const Key& center,
        DISTANCE&& distance_fn = DISTANCE(),
        FILTER&& filter = FILTER()) const {
        return QIteratorKnnHS<Key, T, DISTANCE, FILTER>(
            root_, k, center, std::forward<DISTANCE>(distance_fn), std::forward<FILTER>(filter));
    }

    void check_consistency() {
        // TODO
        assert(true);
    }
};

}  // namespace tinspin
#endif  // TINSPIN_QUADTREE_POINT_H
