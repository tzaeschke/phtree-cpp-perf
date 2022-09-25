/*
 * PHTree.h
 *
 *  Created on: Feb 25, 2016
 *      Author: max
 */

#ifndef SRC_PHTREE_H_
#define SRC_PHTREE_H_

#include <iostream>
#include <vector>
#include <thread>
#include "Entry.h"
#include <thread>

template <unsigned int DIM>
class Node;
template <unsigned int DIM>
class Visitor;
template <unsigned int DIM>
class SizeVisitor;
template <unsigned int DIM, unsigned int WIDTH>
class RangeQueryIterator;
template <unsigned int DIM, unsigned int WIDTH>
class InsertionThreadPool;

template <unsigned int DIM, unsigned int WIDTH>
class PHTree {
	template <unsigned int D>
	friend class SizeVisitor;
	template <unsigned int D, unsigned int W>
	friend std::ostream& operator<<(std::ostream& os, const PHTree<D, W>& tree);
	template <unsigned int D, unsigned int W>
	friend class DynamicNodeOperationsUtil;
	template <unsigned int D, unsigned int W>
	friend class InsertionThreadPool;
public:
	PHTree();
	explicit PHTree(const PHTree<DIM, WIDTH>& other);
	virtual ~PHTree();
	void insert(const Entry<DIM, WIDTH>& e);
	void insert(const std::vector<unsigned long>& values, int id);
	void parallelInsert(const Entry<DIM,WIDTH>& entry);
	void parallelBulkInsert(const std::vector<std::vector<unsigned long>>& values, const std::vector<int>* ids = NULL, size_t nThreads = std::thread::hardware_concurrency());
	void insertHyperRect(const std::vector<unsigned long>& lowerLeftValues, const std::vector<unsigned long>& upperRightValues, int id);
	void bulkInsert(const std::vector<std::vector<unsigned long>>& values, const std::vector<int>& ids);
	void bulkInsert(const std::vector<Entry<DIM,WIDTH>>& entries);

	std::pair<bool,int> lookup(const Entry<DIM, WIDTH>& e) const;
	std::pair<bool,int> lookup(const std::vector<unsigned long>& values) const;
	std::pair<bool,int> lookupHyperRect(const std::vector<unsigned long>& lowerLeftValues, const std::vector<unsigned long>& upperRightValues) const;
	RangeQueryIterator<DIM, WIDTH>* rangeQuery(const Entry<DIM, WIDTH>& lowerLeft, const Entry<DIM, WIDTH>& upperRight) const;
	RangeQueryIterator<DIM, WIDTH>* rangeQuery(const std::vector<unsigned long>& lowerLeftValues, const std::vector<unsigned long>& upperRightValues) const;
	RangeQueryIterator<DIM, WIDTH>* intersectionQuery(const std::vector<unsigned long>& lowerLeftValues, const std::vector<unsigned long>& upperRightValues) const;
	RangeQueryIterator<DIM, WIDTH>* intersectionQuery(const std::vector<unsigned long>& values) const;
	// TODO what exactly to return?
	void parallelIntersectionQuery(const std::vector<std::vector<unsigned long>>& values, size_t nThreads = std::thread::hardware_concurrency()) const;
	RangeQueryIterator<DIM, WIDTH>* inclusionQuery(const std::vector<unsigned long>& lowerLeftValues, const std::vector<unsigned long>& upperRightValues) const;
	RangeQueryIterator<DIM, WIDTH>* inclusionQuery(const std::vector<unsigned long>& values) const;
	void parallelInclusionQuery(const std::vector<std::vector<unsigned long>>& values, size_t nThreads = std::thread::hardware_concurrency()) const;

	void accept(Visitor<DIM>* visitor);

private:
	Node<DIM>* root_;
};

#include <assert.h>
#include "nodes/LHC.h"
#include "util/DynamicNodeOperationsUtil.h"
#include "util/SpatialSelectionOperationsUtil.h"
#include "util/NodeTypeUtil.h"
#include "util/InsertionThreadPool.h"
#include "util/RangeQueryThreadPool.h"

using namespace std;

template <unsigned int DIM, unsigned int WIDTH>
PHTree<DIM, WIDTH>::PHTree() {
	//const unsigned int blocksForFirstSuffix = 1 + ((WIDTH - 1) * DIM - 1) / (8 * sizeof (unsigned long));
	//root_ = NodeTypeUtil<DIM>::template buildNodeWithSuffixes<WIDTH>(0, 1, 1, blocksForFirstSuffix);
    const unsigned int suffixBits = (WIDTH - 1) * DIM; // TODO TZ, this was obviously wrong
    root_ = NodeTypeUtil<DIM>::template buildNodeWithSuffixes<WIDTH>(0, 1, 1, suffixBits);
}

template <unsigned int DIM, unsigned int WIDTH>
PHTree<DIM, WIDTH>::PHTree(const PHTree<DIM, WIDTH>& other) : root_(other.root_) { }

template <unsigned int DIM, unsigned int WIDTH>
PHTree<DIM, WIDTH>::~PHTree() {
	root_->recursiveDelete();
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::insert(const Entry<DIM, WIDTH>& e) {
	#ifdef PRINT
		cout << "inserting: " << e << endl;
	#endif

	DynamicNodeOperationsUtil<DIM, WIDTH>::insert(e, *this);
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::insert(const vector<unsigned long>& values, int id) {
	assert (values.size() == DIM);
	const Entry<DIM, WIDTH> entry(values, id);
	insert(entry);
}


template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::parallelInsert(const Entry<DIM,WIDTH>& entry) {
	DynamicNodeOperationsUtil<DIM,WIDTH>::parallelInsert(entry, this);
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::parallelBulkInsert(const std::vector<std::vector<unsigned long>>& values, const std::vector<int>* ids, size_t nThreads) {
	assert (nThreads > 0);
	InsertionThreadPool<DIM,WIDTH>* pool = new InsertionThreadPool<DIM,WIDTH>(nThreads - 1, values, ids, this);
	pool->joinPool();
	delete pool;
}


template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::bulkInsert(
		const vector<vector<unsigned long>>& values,
		const vector<int>& ids) {
	assert (values.size() == ids.size());


	vector<Entry<DIM,WIDTH>>* entries = new vector<Entry<DIM,WIDTH>>();
	const size_t size = values.size();
	entries->reserve(size);
	for (unsigned i = 0; i < size; ++i) {
		entries->emplace_back(values[i], ids[i]);
	}

	bulkInsert(*entries);
	delete entries;
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::bulkInsert(const vector<Entry<DIM,WIDTH>>& entries) {
	DynamicNodeOperationsUtil<DIM, WIDTH>::bulkInsert(entries, *this);
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::insertHyperRect(
		const vector<unsigned long>& lowerLeftValues,
		const vector<unsigned long>& upperRightValues, int id) {
	assert (DIM % 2 == 0);
	assert (lowerLeftValues.size() == upperRightValues.size());
	assert (lowerLeftValues.size() + upperRightValues.size() == DIM);

	vector<unsigned long> combinedValues(DIM);
	for (unsigned i = 0; i < DIM / 2; ++i) {
		assert (lowerLeftValues[i] <= upperRightValues[i]);
		combinedValues[i] = lowerLeftValues[i];
		combinedValues[i + DIM / 2] = upperRightValues[i];
	}

	insert(combinedValues, id);
}

template <unsigned int DIM, unsigned int WIDTH>
pair<bool,int> PHTree<DIM, WIDTH>:: lookup(const Entry<DIM, WIDTH>& e) const {
	#ifdef PRINT
		cout << "searching: " << e << endl;
	#endif
	return SpatialSelectionOperationsUtil<DIM, WIDTH>::lookup(e, root_, NULL);
}

template <unsigned int DIM, unsigned int WIDTH>
pair<bool,int> PHTree<DIM, WIDTH>::lookup(const std::vector<unsigned long>& values) const {
	const Entry<DIM, WIDTH> entry(values, 0);
	return lookup(entry);
}

template<unsigned int DIM, unsigned int WIDTH>
pair<bool, int> PHTree<DIM, WIDTH>::lookupHyperRect(
		const std::vector<unsigned long>& lowerLeftValues,
		const std::vector<unsigned long>& upperRightValues) const {
	assert (DIM % 2 == 0);
	assert(lowerLeftValues.size() == upperRightValues.size());
	assert(lowerLeftValues.size() + upperRightValues.size() == DIM);

	vector<unsigned long> combinedValues(DIM);
	for (unsigned i = 0; i < DIM / 2; ++i) {
		combinedValues[i] = lowerLeftValues[i];
		combinedValues[i + DIM / 2] = upperRightValues[i];
	}

	return lookup(combinedValues);
}


template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* PHTree<DIM, WIDTH>::rangeQuery(const Entry<DIM, WIDTH>& lowerLeft,
		const Entry<DIM, WIDTH>& upperRight) const {
	// TODO check if lower left and upper right corners are correctly set
	vector<pair<unsigned long, const Node<DIM>*>>* visitedNodes = new vector<pair<unsigned long, const Node<DIM>*>>();
	SpatialSelectionOperationsUtil<DIM, WIDTH>::lookup(lowerLeft, root_, visitedNodes);
	RangeQueryIterator<DIM, WIDTH>* it = new RangeQueryIterator<DIM, WIDTH>(visitedNodes, lowerLeft, upperRight);
	delete visitedNodes;
	return it;
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* PHTree<DIM, WIDTH>::rangeQuery(
		const vector<unsigned long>& lowerLeftValues,
		const vector<unsigned long>& upperRightValues) const {
	const Entry<DIM, WIDTH> lowerLeft(lowerLeftValues, 0);
	const Entry<DIM, WIDTH> upperRight(upperRightValues, 0);
	return rangeQuery(lowerLeft, upperRight);
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* PHTree<DIM, WIDTH>::inclusionQuery(
		const std::vector<unsigned long>& lowerLeftValues,
		const std::vector<unsigned long>& upperRightValues) const {
	assert (DIM % 2 == 0);
	assert ((2 * lowerLeftValues.size() == DIM) && (2 * upperRightValues.size() == DIM));

	vector<unsigned long> lowerLeftHyperRect(DIM);
	vector<unsigned long> upperRightHyperRect(DIM);
	for (unsigned k = 0; k < DIM; ++k) {
		lowerLeftHyperRect[k] = lowerLeftValues[k % (DIM / 2)];
		upperRightHyperRect[k] = upperRightValues[k % (DIM / 2)];
	}

	return rangeQuery(lowerLeftHyperRect, upperRightHyperRect);
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* PHTree<DIM, WIDTH>::inclusionQuery(
		const std::vector<unsigned long>& values) const {
	assert (values.size() == DIM && (DIM % 2 == 0));
	const vector<unsigned long> lowerLeftValues(values.begin(), values.begin() + DIM / 2);
	const vector<unsigned long> upperRightValues(values.begin() + DIM / 2, values.end());
	return inclusionQuery(lowerLeftValues, upperRightValues);
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* PHTree<DIM, WIDTH>::intersectionQuery(
		const vector<unsigned long>& lowerLeftValues,
		const vector<unsigned long>& upperRightValues) const {
	assert (DIM % 2 == 0);
	assert ((2 * lowerLeftValues.size() == DIM) && (2 * upperRightValues.size() == DIM));

	// range query for k-dim hyper rectangles as 2k-dim points:
	// <--- k --->  <-- k -->   <-- k --> <--- k --->
	// (-inf, -inf, ll1, ll2) - (ur1, ur2, +inf, +inf)
	vector<unsigned long> lowerLeftHyperRect(DIM);
	vector<unsigned long> upperRightHyperRect(DIM);
	// set lower half of the values
	for (unsigned k = 0; k <  DIM / 2; ++k) {
		assert (lowerLeftValues[k] <= upperRightValues[k]);
		// with unsigned values 0 is the lowest possible value
		lowerLeftHyperRect[k] = 0;
		upperRightHyperRect[k] = upperRightValues[k];
	}

	// with unsigned values -1 is equal to the highest possible value
	const unsigned long max = (WIDTH == 8 * sizeof (unsigned long))? -1 : (1uL << WIDTH) - 1;
	// set upper half of the values
	for (unsigned k = DIM / 2; k < DIM; ++k) {
		lowerLeftHyperRect[k] = lowerLeftValues[k - DIM / 2];
		upperRightHyperRect[k] = max;
	}

	return rangeQuery(lowerLeftHyperRect, upperRightHyperRect);
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* PHTree<DIM, WIDTH>::intersectionQuery(
		const vector<unsigned long>& values) const {
	assert (values.size() == DIM && (DIM % 2 == 0));
	const vector<unsigned long> lowerLeftValues(values.begin(), values.begin() + DIM / 2);
	const vector<unsigned long> upperRightValues(values.begin() + DIM / 2, values.end());
	return intersectionQuery(lowerLeftValues, upperRightValues);
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::parallelIntersectionQuery(const std::vector<std::vector<unsigned long>>& values, size_t nThreads) const {
	RangeQueryThreadPool<DIM, WIDTH>* pool = new RangeQueryThreadPool<DIM, WIDTH>(nThreads - 1, values, this, intersection_query);
	pool->joinPool();
	// TODO values extracted but what to return?
	delete pool;
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::parallelInclusionQuery(const std::vector<std::vector<unsigned long>>& values, size_t nThreads) const {
	RangeQueryThreadPool<DIM, WIDTH>* pool = new RangeQueryThreadPool<DIM, WIDTH>(nThreads - 1, values, this, inclusion_query);
	pool->joinPool();
	// TODO values extracted but what to return?
	delete pool;
}

template <unsigned int DIM, unsigned int WIDTH>
void PHTree<DIM, WIDTH>::accept(Visitor<DIM>* visitor) {
	(*visitor).template visit<WIDTH>(this);
	root_->accept(visitor, 0, 0);
}

template <unsigned int D, unsigned int W>
ostream& operator <<(ostream& os, const PHTree<D, W> &tree) {
	os << "PH-Tree (dim=" << D << ", value length=" << W << ")" << endl;
	tree.root_->output(os, 0, 0, W);
	return os;
}

#endif /* SRC_PHTREE_H_ */
