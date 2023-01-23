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
    for (int d = 0; d < min.length; d++) {
        if (point[d] < min[d] || point[d] > max[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
static bool isPointEnclosed(const Key& point,
                               const Key& center, double radius) {
    for (int d = 0; d < center.length; d++) {
        if (point[d] < center[d]-radius || point[d] > center[d]+radius) {
            return false;
        }
    }
    return true;
}

template<typename Key>
static bool isPointEqual(const Key& p1, const Key& p2) {
    for (int d = 0; d < p1.length; d++) {
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
    for (int d = 0; d < min.length; d++) {
        if (max[d] < min2[d] || min[d] > max2[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool overlap(const Key& min, const Key& max, const Key& center, double radius) {
    for (int d = 0; d < min.length; d++) {
        if (max[d] < center[d]-radius || min[d] > center[d]+radius) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool overlap(const Key& center, double radius,
                       const Key& center2, double radius2) {
    for (int d = 0; d < center.length; d++) {
        if (center[d]+radius < center2[d]-radius2 || center[d]-radius > center2[d]+radius2) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool isRectEnclosed(const Key& minEnclosed, const Key& maxEnclosed,
                              const Key& minOuter, const Key& maxOuter) {
    for (int d = 0; d < minOuter.length; d++) {
        if (maxOuter[d] < maxEnclosed[d] || minOuter[d] > minEnclosed[d]) {
            return false;
        }
    }
    return true;
}

template<typename Key>
 static bool isRectEnclosed(const Key& minEnclosed, const Key& maxEnclosed,
                              const Key& centerOuter, double radiusOuter) {
    for (int d = 0; d < centerOuter.length; d++) {
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
    for (int d = 0; d < centerOuter.length; d++) {
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
    for (int i = 0; i < p1.length; i++) {
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
    for (int i = 0; i < p.length; i++) {
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
    for (int i = 0; i < center.length; i++) {
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
    for (int i = 0; i < point.length; i++) {
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

template<typename Key, typename T>
class QEntry {
  const Key& point_;
  T value_;

  public: QEntry(const Key& key, T value) {
        point_ = key;
        value_ = value;
    }

    const Key& point() {
        return point;
    }

   const T& value() {
        return value_;
    }

  bool enclosedBy(const Key& min, const Key& max) {
        return isPointEnclosed(point_, min, max);
    }

  bool enclosedBy(const Key& center, double radius) {
        return isPointEnclosed(point_, center, radius);
    }

   bool isExact(QEntry<Key, T> e) {
        return isPointEqual(point_, e.point());
    }


//  String toString() {
//        return "p=" + Arrays.toString(point) + "  v=" + value + " " +
//            System.identityHashCode(this);
//    }

  void setKey(const Key& newPoint) {
        point_ = newPoint;
    }

};

template <typename Key, typename T>
class QEntryDist :public QEntry<Key, T> { // TODO inheritance???
  double distance_;

public: QEntryDist(QEntry<Key, T> e, double dist) :

      QEntry<Key, T>(e.point(), e.value()) {
        distance_ = dist;
    }

   double dist() {
        return distance_;
    }

//  static final QEntryComparator COMP = new QEntryComparator();
//
//    static class QEntryComparator implements Comparator<QEntryDist<?>> {
//
//        /**
//	    * Compares the two specified MBRs according to
//	    * the sorting dimension and the sorting co-ordinate for the dimension
//	     * of this Comparator.
//	    *
//	    * @param o1 the first SpatialPoint
//	    * @param o2 the second SpatialPoint
//	    * @return a negative integer, zero, or a positive integer as the
//	    *         first argument is less than, equal to, or greater than the
//	    *         second.
//         */
//        @Override
//	    public int compare(QEntryDist<?> o1, QEntryDist<?> o2) {
//            double d = o1.dist() - o2.dist();
//            return d < 0 ? -1 : (d > 0 ? 1 : 0);
//        }
//    }

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
        subs_.add(subNode);
    }

    QNode<Key, T>* tryPut(QEntry<Key, T>& e, size_t maxNodeSize, bool enforceLeaf) {
        if (DEBUG && !e.enclosedBy(center_, radius_)) {
            std::cerr << "ERROR: e=" + e.point() << " center/radius=" << center_ << "/" << radius_ << std::endl;
        }

        // traverse subs?
        if (values_ == nullptr) {
            return getOrCreateSub(e);
        }

        // add if:
        // a) we have space
        // b) we have maxDepth
        // c) elements are equal (work only for n=1, avoids splitting
        //    in cases where splitting won't help. For n>1 the
        //    local limit is (temporarily) violated.
        if (values_.size() < maxNodeSize || enforceLeaf || e.isExact(values_.get(0))) {
            values_.add(e);
            return nullptr;
        }

        // split
        std::vector<QEntry<Key, T>> vals = std::move(values_);  // TODO avoid move and erase later?
        values_.clear(); // = nullptr;
        values_.shrink_to_fit();
        assert(subs_.empty());
        for (int i = 0; i < vals.size(); i++) {
            QEntry<Key, T> e2 = vals.get(i);
            QNode<Key, T> sub = getOrCreateSub(e2);
            while (sub != nullptr) {
                // This may recurse if all entries fall
                // into the same subnode
                sub = (QNode<Key, T>)sub.tryPut(e2, maxNodeSize, false);
            }
        }
        return getOrCreateSub(e);
    }

  private:
    QNode<Key, T>* getOrCreateSub(const QEntry<Key, T>& e) {
        QNode<Key, T>* n = findSubNode(e.point());
        if (n == nullptr) {
            n = createSubForEntry(e.point());
            subs_.add(n);
        }
        return n;
    }

    QNode<Key, T>* createSubForEntry(const Key& p) {
        Key centerSub = new double[center_.length];
        // This ensures that the subsnodes completely cover the area of
        // the parent node.
        double radiusSub = radius_ / 2.0;
        for (int d = 0; d < center_.length; d++) {
            if (p[d] >= center_[d]) {
                centerSub[d] = center_[d] + radiusSub;
            } else {
                centerSub[d] = center_[d] - radiusSub;
            }
        }
        return new QNode<Key, T>(centerSub, radiusSub);
    }

    /**
	 * The subnode position has reverse ordering of the point's
	 * dimension ordering. Dimension 0 of a point is the highest
	 * ordered bit in the position.
	 * @param p point
	 * @return subnode position
     */
    QNode<Key, T>* findSubNode(const Key& p) {
        for (int i = 0; i < subs_.size(); i++) {
            QNode<Key, T>* n = subs_.get(i);
            if (isPointEnclosed(p, n->center_, n->radius_)) {
                return n;
            }
        }
        return nullptr;
    }

  public:
    QEntry<Key, T> remove(QNode<Key, T>* parent, const Key& key, int maxNodeSize) {
        if (values_ == nullptr) {
            QNode<Key, T> sub = findSubNode(key);
            if (sub != nullptr) {
                return sub.remove(this, key, maxNodeSize);
            }
            return nullptr;
        }

        for (int i = 0; i < values_.size(); i++) {
            QEntry<Key, T> e = values_.get(i);
            if (isPointEqual(e.point(), key)) {
                values_.remove(i);
                // TODO provide threshold for re-insert
                // i.e. do not always merge.
                if (parent != nullptr) {
                    parent->checkAndMergeLeafNodes(maxNodeSize);
                }
                return e;
            }
        }
        return nullptr;
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
            QEntry<Key, T> ret = sub->update(
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

        for (int i = 0; i < values_.size(); i++) {
            QEntry<Key, T> e = values_.get(i);
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
        for (int i = 0; i < subs_.size(); i++) {
            if (subs_.get(i).values_ == nullptr) {
                // can't merge directory nodes.
                return;
            }
            nTotal += subs_.get(i).values_.size();
            if (nTotal > maxNodeSize) {
                // too many children
                return;
            }
        }

        // okay, let's merge
        assert(values_.empty());
        values_.reserve(nTotal);
        //values_ = new ArrayList<>(nTotal);
        for (int i = 0; i < subs_.size(); i++) {
            values_.addAll(subs_.get(i).values_);
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

    QEntry<Key, T> getExact(const Key& key) {
        if (values_ == nullptr) {
            QNode<Key, T> sub = findSubNode(key);
            if (sub != nullptr) {
                return sub.getExact(key);
            }
            return nullptr;
        }

        for (int i = 0; i < values_.size(); i++) {
            QEntry<Key, T> e = values_.get(i);
            if (isPointEqual(e.point(), key)) {
                return e;
            }
        }
        return nullptr;
    }

    auto& getEntries() {
        return values_;
    }

    auto getChildIterator() {
        if (values_ != nullptr) {
            return values_.iterator();
        }
        return subs_.iterator();
    }

//    public: String toString() {
//        return "center/radius=" + Arrays.toString(center) + "/" + radius + " " +
//            System.identityHashCode(this);
//    }

    void checkNode(QStats s, QNode<Key, T> parent, int depth) {
        if (depth > s.maxDepth) {
            s.maxDepth = depth;
        }
        s.nNodes++;

        if (parent != nullptr) {
            if (!isRectEnclosed(
                    center_, radius_, parent.center_, parent.radius * EPS_MUL)) {
                for (int d = 0; d < center_.length; d++) {
                    //					if ((centerOuter[d]+radiusOuter) / (centerEnclosed[d]+radiusEnclosed) < 0.9999999 || 							(centerOuter[d]-radiusOuter) / (centerEnclosed[d]-radiusEnclosed) > 1.0000001) { 						return false;
                    //					}
                    std::cout <<
                        "Outer: " <<  parent.radius_ <<  " " <<  parent.center_ << std::endl;
                    std::cout << "Child: " <<  radius_ <<  " " <<  center_ << std::endl;
                    std::cout <<
                        (parent.center_[d] +  parent.radius_) <<  " vs " <<  (center_[d] +  radius_) << std::endl;
                    std::cout <<
                        "r=" <<  (parent.center_[d] +  parent.radius_) / (center_[d] +  radius_) << std::endl;
                    std::cout <<
                        (parent.center_[d] - parent.radius_) <<  " vs " <<  (center_[d] - radius_) << std::endl;
                    std::cout <<
                        "r=" <<  (parent.center_[d] - parent.radius_) / (center_[d] - radius_) << std::endl;
                }
                assert(false);
            }
        }
        if (values_ != nullptr) {
            for (int i = 0; i < values_.size(); i++) {
                QEntry<Key, T> e = values_.get(i);
                if (!isPointEnclosed(e.point(), center_, radius_ * EPS_MUL)) {
                    std::cout << "Node: " <<  radius_ <<  " " <<  center_ << std::endl;
                    std::cout << "Child: " <<  e.point() << std::endl;
                    for (int d = 0; d < center_.length; d++) {
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
            for (int i = 0; i < subs_.size(); i++) {
                QNode<Key, T> n = subs_.get(i);
                n.checkNode(s, this, depth + 1);
            }
        }
    }

    bool isLeaf() {
        return values_ != nullptr;
    }

    auto& getChildNodes() {
        return subs_;
    }
};
}

 template<typename T>
class QuadTree {
    static const int MAX_DEPTH = 50;
    using Key = PhPointD<3>;

    static const int DEFAULT_MAX_NODE_SIZE = 10;

    const size_t dims;
    const size_t maxNodeSize;
    QNode<Key, T> root = nullptr;
    size_t size_ = 0;

  public:
    QuadTree(size_t dims, size_t maxNodeSize = DEFAULT_MAX_NODE_SIZE) {
        if (DEBUG) {
            std::cout << "Warning: DEBUG enabled" << std::endl; // TODO
        }
        this->dims = dims;
        this->maxNodeSize = maxNodeSize;
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
    void insert(const Key& key, const T& value) {
        size_++;
        QEntry<Key, T> e(key, value); // TODO std::move into node
        if (root == nullptr) {
            initializeRoot(key);
        }
        ensureCoverage(e);
        auto* r = root;
        int depth = 0;
        while (r != nullptr) {  // TODO backport, r is always a QNode
            r = r->tryPut(e, maxNodeSize, depth++ > MAX_DEPTH);
        }
    }

  private:
    void initializeRoot(const Key& key) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (int d = 0; d < dims; d++) {
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
        for (int d = 0; d < dims; d++) {
            center[d] = key[d] > 0 ? maxDistOrigin : -maxDistOrigin;
            //			max[d] = key[d] < 0 ? 0 : (maxDistOrigin*2);
        }
        root = new QNode<Key, T>(center, maxDistOrigin);
    }

    /**
	 * Check whether a given key exists.
	 * @param key the key to check
	 * @return true iff the key exists
     */
  public:
    bool containsExact(const Key& key) {
        if (root == nullptr) {
            return false;
        }
        return root.getExact(key) != nullptr;
    }

    /**
	 * Get the value associates with the key.
	 * @param key the key to look up
	 * @return the value for the key or 'nullptr' if the key was not found
     */
    public: T queryExact(const Key& key) {
        if (root == nullptr) {
            return nullptr;
        }
        QEntry<Key, T> e = root.getExact(key);
        return e == nullptr ? nullptr : e.value();
    }

    /**
	 * Remove a key.
	 * @param key key to remove
	 * @return the value associated with the key or 'nullptr' if the key was not found
     */
    public: T remove(const Key& key) {
        if (root == nullptr) {
            if (DEBUG) {
                std::cerr <<"Failed remove 1: " << key << std::endl;
            }
            return nullptr;
        }
        QEntry<Key, T> e = root.remove(nullptr, key, maxNodeSize);
        if (e == nullptr) {
            if (DEBUG) {
                std::cerr << "Failed remove 2: " << key << std::endl;
            }
            return nullptr;
        }
        size_--;
        return e.value();
    }

    /**
	 * Reinsert the key.
	 * @param oldKey old key
	 * @param newKey new key
	 * @return the value associated with the key or 'nullptr' if the key was not found.
     */
  public: T& update(const Key& oldKey, const Key& newKey) {
        if (root == nullptr) {
            return nullptr;
        }
        bool requiresReinsert = false;
        QEntry<Key, T> e =
            root.update(nullptr, oldKey, newKey, maxNodeSize, requiresReinsert, 0, MAX_DEPTH);
        if (e == nullptr) {
            // not found
            if (DEBUG) {
                std::cout << "Failed reinsert 1: " << oldKey << "/" << newKey << std::endl;
            }
            return nullptr;
        }
        if (requiresReinsert) {
            if (DEBUG) {
                std::cout << "Failed reinsert 2: " << oldKey << "/" << newKey << std::endl;
            }
            // does not fit in root node...
            ensureCoverage(e);
            auto* r = root;
            int depth = 0;
            while (r != nullptr) {  // TODO backport, r is always a QNode
                r = r->tryPut(e, maxNodeSize, depth++ > MAX_DEPTH);
            }
        }
        return e.value();
    }

    /**
	 * Ensure that the tree covers the entry.
	 * @param e Entry to cover.
     */
    private: void ensureCoverage(QEntry<Key, T> e) {
        const Key& p = e.point();
        while (!e.enclosedBy(root.getCenter(), root.getRadius())) {
            const Key& center = root.getCenter();
            double radius = root.getRadius();
            Key center2{};
            double radius2 = radius * 2;
            // TODO use DIM?
            for (int d = 0; d < center.size(); d++) {
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
            root = new QNode<Key, T>(center2, radius2, root);
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
    public: void clear() {
        size_ = 0;
        root = nullptr;
    }

    /**
	 * Query the tree, returning all points in the axis-aligned rectangle between 'min' and 'max'.
	 * @param min lower left corner of query
	 * @param max upper right corner of query
	 * @return all entries in the rectangle
     */
    public: auto query(const Key& min, const Key& max) {
        return QIterator<Key, T>(this, min, max);
    }

    /**
	 * Resettable query iterator.
	 *
	 * @param <T> Value type
     */
  public:
    template <typename Key, typename T>
    class QIterator { //implements QueryIterator<PointEntry<T>> {
        using IterT = decltype(std::vector<QNode<Key, T>>().begin());
        const QuadTree<T> tree;
        std::stack <IterT> stack; // TODO backport, Why is this a Deque????????
        QEntry<Key, T> next = nullptr;
        Key min;
        Key max;

        QIterator(QuadTree<T> tree, const Key& min, const Key& max) : tree{tree}, stack{}, next{nullptr}, min(min), max(max){
            if (tree.root != nullptr) {
                stack.push(tree.root.getChildIterator());
                findNext();
            }
        }

         private: void findNext() {
            while (!stack.isEmpty()) {
                auto& it = stack.peek();
                while (it.hasNext()) {
                    Object o = it.next();
                    if (o instanceof QNode) {
                        QNode<Key, T>* node = (QNode<Key, T>)o;
                        if (overlap(min, max, node->getCenter(), node->getRadius())) {
                            it = node->getChildIterator();
                            stack.push(it);
                        }
                        continue;
                    }
                    QEntry<Key, T> e = (QEntry<Key, T>)o;
                    if (e.enclosedBy(min, max)) {
                        next = e;
                        return;
                    }
                }
                stack.pop();
            }
            next = nullptr;
        }

        public: bool hasNext() {
            return next != nullptr;
        }

        public: const QEntry<Key, T>& next2() {
            assert(hasNext());
            QEntry<Key, T> ret = next;
            findNext();
            return ret;
        }
    };

    public: std::vector<QEntryDist<Key, T>>
    knnQuery(const Key& center, int k) {
        if (root == nullptr) {
            return std::vector<QEntryDist<Key, T>>{};
        }
        Comparator<QEntry<Key, T>> comp = (QEntry<Key, T> point1, QEntry<Key, T> point2)->{
            double deltaDist =
                distance(center, point1.point()) - distance(center, point2.point());
            return deltaDist < 0 ? -1 : (deltaDist > 0 ? 1 : 0);
        };
        double distEstimate = distanceEstimate(root, center, k, comp);
        std::vector<QEntryDist<Key, T>> candidates{};
        candidates.reserve(k);
        while (candidates.size() < k) {
            candidates.clear();
            rangeSearchKNN(root, center, candidates, k, distEstimate);
            distEstimate *= 2;
        }
        return candidates;
    }

    private:
      template <typename COMP>
      double distanceEstimate(QNode<Key, T>* node, const Key& point, size_t k, const COMP& comp) {
        if (node->isLeaf()) {
            // This is a leaf that would contain the point.
            int n = node->getEntries().size();
            // Create a copy!
            std::vector<QEntry<Key, T>> data(node->getEntries());
            std::sort(data.begin(), data.end(), comp);
            int pos = n < k ? n : k;
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
            for (int i = 0; i < nodes.size(); i++) {
                QNode<Key, T> sub = nodes.get(i);
                if (isPointEnclosed(point, sub.getCenter(), sub.getRadius())) {
                    return distanceEstimate(sub, point, k, comp);
                }
            }
            // okay, this directory node contains the point, but none of the leaves does.
            // We just return the size of this node, because all it's leaf nodes should
            // contain more than enough candidate in proximity of 'point'.
            return node->getRadius() * Math.sqrt(point.length);
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
            for (int i = 0; i < points.size(); i++) {
                QEntry<Key, T> p = points.get(i);
                double dist = distance(center, p.point());
                if (dist < maxRange) {
                    candidates.add(new QEntryDist<Key, T>(p, dist));
                }
            }
            maxRange = adjustRegionKNN(candidates, k, maxRange);
        } else {
            std::vector<QNode<Key, T>>& nodes = node->getChildNodes();
            for (int i = 0; i < nodes.size(); i++) {
                QNode<Key, T> sub = nodes.get(i);
                if (sub != nullptr &&
                    distToRectNode(center, sub.getCenter(), sub.getRadius()) < maxRange) {
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

        // use stored distances instead of recalcualting them
        candidates.sort(QEntryDist.COMP);
        while (candidates.size() > k) {
            candidates.remove(candidates.size() - 1);
        }

        double range = candidates.get(candidates.size() - 1).dist();
        return range;
    }

  private:
    class QQueryIteratorKNN { //implements QueryIteratorKNN<PointEntryDist<T>> {
        using IterT = decltype(std::vector<QEntryDist<Key, T>>().begin());
        IterT it;

      public:
        QQueryIteratorKNN(const Key& center, int k) {
            it = knnQuery(center, k).begin();
        }

        public: bool hasNext() {
            return it.hasNext();
        }

        public: const QEntryDist<Key, T>& next() {
            return it.next();
        }
    };

    /**
	 * Returns a printable list of the tree.
	 * @return the tree as String
     */
//    public: String
//    toStringTree() {
//        StringBuilder sb = new StringBuilder();
//        if (root == nullptr) {
//            sb.append("empty tree");
//        } else {
//            toStringTree(sb, root, 0, 0);
//        }
//        return sb.toString();
//    }
//
//     private: void toStringTree(
//        StringBuilder sb, QNode<Key, T> node, int depth, int posInParent) {
//        Iterator < ? > it = node->getChildIterator();
//        String prefix = "";
//        for (int i = 0; i < depth; i++) {
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
//                QNode<Key, T> sub = (QNode<Key, T>)o;
//                toStringTree(sb, sub, depth + 1, pos);
//            } else if (o instanceof QEntry) {
//                QEntry<Key, T> e = (QEntry<Key, T>)o;
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
//            (root == nullptr ? "nullptr" : (Arrays.toString(root.getCenter()) + "/" + root.getRadius()));
//    }

    public: QStats getStats() {
        QStats s{};
        if (root != nullptr) {
            root.checkNode(s, nullptr, 0);
        }
        return s;
    }

    public: int
    getDims() {
        return dims;
    }

    public:
      auto begin() {
        if (root == nullptr) {
            return query(new double[dims], new double[dims]);
        }
        // return query(root.);
        // TODO
        assert(false);
    }

    public: QQueryIteratorKNN queryKNN(const Key& center, int k) {
        return new QQueryIteratorKNN(center, k);
    }

    public: int getNodeCount() {
        return getStats().getNodeCount();
    }

    public: int getDepth() {
        return getStats().getMaxDepth();
    }
};

}
#endif  // TINSPIN_QUADTREE_POINT_H
