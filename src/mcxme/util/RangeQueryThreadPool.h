/*
 * RangeQueryThreadPool.h
 *
 *  Created on: Aug 19, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_RANGEQUERYTHREADPOOL_H_
#define SRC_UTIL_RANGEQUERYTHREADPOOL_H_


#include <thread>
#include <queue>
#include <functional>
#include <vector>
#include <atomic>
#include "../Entry.h"
#include "../util/ResultStorage.h"

template <unsigned int DIM, unsigned int WIDTH>
class PHTree;

template <unsigned int DIM, unsigned int WIDTH>
class ResultStorage;

enum QueryType {
	intersection_query,
	inclusion_query
};

template <unsigned int DIM, unsigned int WIDTH>
class RangeQueryThreadPool {
public:

	RangeQueryThreadPool(size_t nAdditionalThreads,
			const std::vector<std::vector<unsigned long>>& ranges,
			const PHTree<DIM, WIDTH>* tree, QueryType type);
	~RangeQueryThreadPool();
	void joinPool();

private:
	QueryType type_;
	size_t nThreads_;
	std::vector<std::thread> threads_;
	std::vector<ResultStorage<DIM, WIDTH>> startStorage_;
	const std::vector<std::vector<unsigned long>>& ranges_;
	const PHTree<DIM, WIDTH>* tree_;

	void processNext(size_t threadIndex);
	RangeQueryIterator<DIM, WIDTH>* getIteratorByType(size_t index);
	void clearResults();
};


using namespace std;

#include <stdexcept>
#include "../PHTree.h"

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryThreadPool<DIM, WIDTH>::RangeQueryThreadPool(size_t nAdditionalThreads,
			const std::vector<std::vector<unsigned long>>& ranges,
			const PHTree<DIM, WIDTH>* tree, QueryType type) :
			type_(type), nThreads_(nAdditionalThreads + 1),
			threads_(), startStorage_(nAdditionalThreads + 1), ranges_(ranges), tree_(tree) {
	threads_.reserve(nAdditionalThreads);
	for (unsigned tCount = 0; tCount < nAdditionalThreads; ++tCount) {
		threads_.emplace_back(&RangeQueryThreadPool<DIM,WIDTH>::processNext, this, tCount);
	}
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryThreadPool<DIM, WIDTH>::~RangeQueryThreadPool() {
	for (auto &t : threads_) {
		t.join();
	}

	clearResults();
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryThreadPool<DIM, WIDTH>::joinPool() {
	processNext(nThreads_ - 1);
}


template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>* RangeQueryThreadPool<DIM, WIDTH>::getIteratorByType(size_t index) {
	switch (type_) {
	case intersection_query: return tree_->intersectionQuery(ranges_[index]);
	case inclusion_query: return tree_->inclusionQuery(ranges_[index]);
	default: throw runtime_error("unknown query type");
	}
}


template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryThreadPool<DIM, WIDTH>::clearResults() {
	for (auto &start : startStorage_) {
		ResultStorage<DIM, WIDTH>* storage = start.nextStorage_;
		while (storage != NULL) {
			ResultStorage<DIM, WIDTH>* nextStorage = storage->nextStorage_;
			delete storage;
			storage = nextStorage;
		}
	}
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryThreadPool<DIM, WIDTH>::processNext(size_t threadIndex) {
	const size_t chunkSize = 1 + ranges_.size() / nThreads_;
	const size_t start = chunkSize * threadIndex;
	const size_t end = min(chunkSize * (threadIndex + 1), ranges_.size());
	ResultStorage<DIM, WIDTH>* storage = &startStorage_[threadIndex];

	for (size_t i = start; i < end; ++i) {
		RangeQueryIterator<DIM, WIDTH>* it = getIteratorByType(i);
		while (it->hasNext()) {
			const Entry<DIM, WIDTH> entry = it->next();
			storage = storage->add(entry);
		}

		delete it;
	}
}

#endif /* SRC_UTIL_RANGEQUERYTHREADPOOL_H_ */
