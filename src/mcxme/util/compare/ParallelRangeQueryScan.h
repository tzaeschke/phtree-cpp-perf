/*
 * ParallelRangeQueryScan.h
 *
 *  Created on: Aug 21, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_PARALLELRANGEQUERYSCAN_H_
#define SRC_UTIL_PARALLELRANGEQUERYSCAN_H_


#include <thread>
#include <vector>
#include <mutex>

class ParallelRangeQueryScan {
public:

	ParallelRangeQueryScan(size_t nAdditionalThreads,
			const std::vector<std::vector<unsigned long>>& entries,
			const std::vector<unsigned long>& range);
	~ParallelRangeQueryScan();
	void joinWork();
	void finishWork();

	const std::vector<std::vector<unsigned long>>* getResult() const;

private:
	size_t nThreads_;
	std::mutex m_;
	std::vector<std::thread> threads_;
	const std::vector<unsigned long>& range_;
	const std::vector<std::vector<unsigned long>>& entries_;
	std::vector<std::vector<unsigned long>> result_;
	void process(size_t threadIndex);
};


using namespace std;

#include <assert.h>

ParallelRangeQueryScan::ParallelRangeQueryScan(size_t nAdditionalThreads,
			const vector<vector<unsigned long>>& entries,
			const vector<unsigned long>& range) :
			nThreads_(nAdditionalThreads + 1), m_(), threads_(), range_(range),
			entries_(entries), result_() {

	threads_.reserve(nAdditionalThreads);
	for (unsigned tCount = 0; tCount < nAdditionalThreads; ++tCount) {
		threads_.emplace_back(&ParallelRangeQueryScan::process, this, tCount);
	}
}

ParallelRangeQueryScan::~ParallelRangeQueryScan() { }

const std::vector<std::vector<unsigned long>>* ParallelRangeQueryScan::getResult() const {
	return &result_;
}

void ParallelRangeQueryScan::joinWork() {
	process(threads_.size());
}

void ParallelRangeQueryScan::finishWork() {
	for (auto &t : threads_) {
		t.join();
	}
}


void ParallelRangeQueryScan::process(size_t threadIndex) {
	assert (entries_.size() > 0);

	const size_t chunkSize = 1 + entries_.size() / nThreads_;
	const size_t start = chunkSize * threadIndex;
	const size_t end = min(chunkSize * (threadIndex + 1), entries_.size());
	const size_t dim = entries_[0].size();
	assert (dim > 0);
	assert (range_.size() == dim);

	for (unsigned i = start; i < end; ++i) {
		assert (entries_[i].size() == dim);
		// check if the current entry is within the range
		bool valid = true;
		for (unsigned j = 0; j < dim / 2 && valid; ++j) {
			const unsigned long dimValue = entries_[i][j];
			valid = (0 <= dimValue) && (dimValue <= range_[j + dim / 2]);
		}

		for (unsigned j = dim / 2; j < dim && valid; ++j) {
			const unsigned long dimValue = entries_[i][j];
			valid = (range_[j - dim / 2] <= entries_[i][j]) && (dimValue <= (-1uL));
		}

		if (valid) {
			unique_lock<mutex> lk(m_);
			result_.push_back(entries_[i]);
		}
	}
}

#endif /* SRC_UTIL_PARALLELRANGEQUERYSCAN_H_ */
