/*
 * ThreadPool.h
 *
 *  Created on: Jul 18, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_INSERTIONTHREADPOOL_H_
#define SRC_UTIL_INSERTIONTHREADPOOL_H_

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <vector>
#include <atomic>
#include <boost/thread/barrier.hpp>
#include <boost/thread/shared_mutex.hpp>
#include "../util/DeletedNodes.h"
#include "../util/EntryTreeMap.h"

template <unsigned int DIM, unsigned int WIDTH>
class PHTree;
template <unsigned int DIM, unsigned int WIDTH>
class Entry;
template <unsigned int DIM, unsigned int WIDTH>
class EntryBufferPool;


enum InsertionOrder {
	sequential_entries,
	range_per_thread, // DEFAULT
	sequential_ranges
};

enum InsertionApproach {
	optimistic_locking,
	buffered_bulk // DEFAULT
};

template <unsigned int DIM, unsigned int WIDTH>
class InsertionThreadPool {
public:

	InsertionThreadPool(size_t furtherThreads,
			const std::vector<std::vector<unsigned long>>& values,
			const std::vector<int>* ids, PHTree<DIM, WIDTH>* tree);
	~InsertionThreadPool();
	void joinPool();

	const size_t fixRangeSize = 100;
	static InsertionOrder order_;
	static InsertionApproach approach_;

	static unsigned long nFlushPhases;

private:

	static const size_t INITIAL_JITTER_NANOS_FACTOR = 500;

	bool syncPhaseRequired_;
	std::atomic<unsigned int> i_;
	size_t nThreads_;
	std::atomic<size_t> nRemainingThreads_;
	boost::shared_mutex createBarriersMutex_;
	boost::barrier* poolFlushBarrier_;
	std::vector<std::thread> threads_;
	std::vector<DeletedNodes<DIM>> deletedNodes_;
	std::vector<EntryTreeMap<DIM,WIDTH>> entryMaps_;
	std::vector<std::vector<double>> nanosPerEntryPerThread_;
	const std::vector<std::vector<unsigned long>>& values_;
	const std::vector<int>* ids_;
	PHTree<DIM, WIDTH>* tree_;
	EntryBufferPool<DIM, WIDTH>* pool_;

	void processNext(size_t threadIndex);
	inline double insertBySelectedStrategy(size_t entryIndex, size_t threadIndex);

	inline void handlePoolFlushSync(size_t threadIndex, bool lastFlush);
	inline void handleDoneBySelectedStrategy(size_t threadIndex);

};

#include <iostream>
#include <string>
#include <fstream>

#include "../Entry.h"
#include "../util/DynamicNodeOperationsUtil.h"
#include "../util/NodeTypeUtil.h"
#include "../util/EntryBufferPool.h"

using namespace std;

template <unsigned int DIM, unsigned int WIDTH>
InsertionOrder InsertionThreadPool<DIM, WIDTH>::order_ = range_per_thread;
template <unsigned int DIM, unsigned int WIDTH>
InsertionApproach InsertionThreadPool<DIM, WIDTH>::approach_ = buffered_bulk;
template <unsigned int DIM, unsigned int WIDTH>
unsigned long InsertionThreadPool<DIM, WIDTH>::nFlushPhases = 0;

template <unsigned int DIM, unsigned int WIDTH>
InsertionThreadPool<DIM, WIDTH>::InsertionThreadPool(size_t furtherThreads,
		const vector<vector<unsigned long>>& values,
		const vector<int>* ids, PHTree<DIM, WIDTH>* tree)
		: syncPhaseRequired_(false), i_(0), nThreads_(furtherThreads + 1), createBarriersMutex_(),
		  poolFlushBarrier_(NULL), nanosPerEntryPerThread_(furtherThreads + 1), deletedNodes_(furtherThreads + 1),
		  entryMaps_(furtherThreads + 1), values_(values),
		  ids_(ids), tree_(tree), pool_(NULL) {
	assert (values.size() > 0);
	// create the biggest possible root node so there is no need to synchronize access on the root
	Node<DIM>* oldRoot = tree->root_;
	assert (oldRoot->getNumberOfContents() == 0);
	Node<DIM>* newRoot = NodeTypeUtil<DIM>::copyIntoLargerNode(1uL << DIM, oldRoot);
	tree->root_ = newRoot;
	delete oldRoot;

	poolFlushBarrier_ = new boost::barrier(nThreads_);
	pool_ = new EntryBufferPool<DIM,WIDTH>(); // TODO only create if needed

	DynamicNodeOperationsUtil<DIM,WIDTH>::nThreads = nThreads_;
	nRemainingThreads_ = nThreads_;
	threads_.reserve(furtherThreads);
	for (unsigned tCount = 0; tCount < furtherThreads; ++tCount) {
		threads_.emplace_back(&InsertionThreadPool<DIM,WIDTH>::processNext, this, tCount);
	}
}

template <unsigned int DIM, unsigned int WIDTH>
InsertionThreadPool<DIM, WIDTH>::~InsertionThreadPool() {
	for (auto &t : threads_) {
		t.join();
	}

	assert (nRemainingThreads_ == 0);
	delete poolFlushBarrier_;
	delete pool_;

	// TODO remove:
	/*string path = "./plot/data/timeseries-" + to_string(nThreads_) + ".dat";
	ofstream* dataFile = new ofstream();
	dataFile->open(path.c_str(), ofstream::out | ofstream::trunc);
	size_t nEntries = nanosPerEntryPerThread_[nThreads_-1].size();
	for (unsigned e = 0; e < nEntries; ++e) {
		(*dataFile) << e;
		for (unsigned t = 0; t < nThreads_; ++t) {
			(*dataFile) << "\t" << nanosPerEntryPerThread_[t][e];
		}
		(*dataFile) << endl;
	}
	delete dataFile;*/
	// shrink the root node again
	size_t rootContents = tree_->root_->getNumberOfContents();
	assert (rootContents > 0);
	Node<DIM>* oldRoot = tree_->root_;
	Node<DIM>* newRoot = NodeTypeUtil<DIM>::copyIntoLargerNode(rootContents, oldRoot);
	tree_->root_ = newRoot;
	delete oldRoot;
}

template <unsigned int DIM, unsigned int WIDTH>
void InsertionThreadPool<DIM, WIDTH>::joinPool() {
	processNext(threads_.size());
}

template <unsigned int DIM, unsigned int WIDTH>
void InsertionThreadPool<DIM, WIDTH>::handlePoolFlushSync(size_t threadIndex, bool lastFlush) {

	// gather all threads
	const bool responsibleForState = poolFlushBarrier_->wait();
	if (responsibleForState) {
		pool_->prepareFullDeallocate();
	}

	// wait until the pool is prepared for the flush
	poolFlushBarrier_->wait();

	pool_->doFullDeallocatePart(threadIndex, nThreads_);
	deletedNodes_[threadIndex].deleteAll();
	entryMaps_[threadIndex].clearMap();

	if (responsibleForState) {
		++nFlushPhases;
		syncPhaseRequired_ = false;
	}

	if (lastFlush) { --nRemainingThreads_; }

	// wait until every thread finished flushing its part
	poolFlushBarrier_->wait();
	if (responsibleForState) { // TODO maybe possible to skip this step
		pool_->finishFullDeallocate();
	}

	// wait until the pool is ready to restart
	poolFlushBarrier_->wait();
}

template <unsigned int DIM, unsigned int WIDTH>
void InsertionThreadPool<DIM, WIDTH>::handleDoneBySelectedStrategy(size_t threadIndex) {
	switch (approach_) {
	case buffered_bulk:
		handlePoolFlushSync(threadIndex, true);
		while (nRemainingThreads_ > 0) {
			handlePoolFlushSync(threadIndex, false);
		}
		break;
	case optimistic_locking:
		--nRemainingThreads_;
		break;
	default:
		break;
	}
}

template <unsigned int DIM, unsigned int WIDTH>
double InsertionThreadPool<DIM, WIDTH>::insertBySelectedStrategy(size_t entryIndex, size_t threadIndex) {

	struct timespec start, finish;
	clock_gettime(CLOCK_MONOTONIC, &start);

	const int id = (ids_)? (*ids_)[entryIndex] : entryIndex;
	const Entry<DIM, WIDTH>* entry = entryMaps_[threadIndex].createEntry(values_[entryIndex], id);
	switch (approach_) {
	case optimistic_locking:
		DynamicNodeOperationsUtil<DIM, WIDTH>::parallelInsert(*entry, *tree_);
		break;
	case buffered_bulk:
		bool success = false;
		while (!success) {
			if (deletedNodes_[threadIndex].full()) {syncPhaseRequired_ = true; }
			if (syncPhaseRequired_) { handlePoolFlushSync(threadIndex, false); }
			success = DynamicNodeOperationsUtil<DIM, WIDTH>::parallelBulkInsert(
					*entry, *tree_, *pool_,
					deletedNodes_[threadIndex],
					entryMaps_[threadIndex]);
			if (!success) { syncPhaseRequired_ = true; }
		}
		break;
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	const double elapsedNanos = finish.tv_nsec - start.tv_nsec;
	return elapsedNanos;
}

template <unsigned int DIM, unsigned int WIDTH>
void InsertionThreadPool<DIM, WIDTH>::processNext(size_t threadIndex) {


	const size_t size = values_.size();
	switch (order_) {
	case sequential_entries:
		{
			size_t i = 0;
			while (i < size) {
				i = i_++;
				if (i < size) {
					insertBySelectedStrategy(i, threadIndex);
				}
			}
		}
		break;
	case range_per_thread:
		{
			// TODO misses last entries for size mod thread != 0
			const size_t start = size * threadIndex / nThreads_;
			const size_t end = min(size * (threadIndex + 1) / nThreads_, size);
//			nanosPerEntryPerThread_[threadIndex].resize(end - start);
			for (size_t i = start; i < end; ++i) {
				double nanos = insertBySelectedStrategy(i, threadIndex);
//				nanosPerEntryPerThread_[threadIndex][i - start] = nanos;
			}
		}
		break;
	case sequential_ranges:
		{
			size_t start, end;
			start = 0;
			while (start < size) {
				const size_t blockIndex = i_++;
				start = blockIndex * fixRangeSize;
				end = min(start + fixRangeSize, size);
				for (size_t i = start; i < end; ++i) {
					insertBySelectedStrategy(i, threadIndex);
				}
			}
		}
		break;
	default: throw "unknown order";
	}

	handleDoneBySelectedStrategy(threadIndex);
#ifdef PRINT
		cout << "thread (ID: " << threadIndex << ") finished" << endl;
#endif
}


#endif /* SRC_UTIL_INSERTIONTHREADPOOL_H_ */
