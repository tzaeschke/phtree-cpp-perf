/*
 * ResultStorage.h
 *
 *  Created on: Aug 20, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_RESULTSTORAGE_H_
#define SRC_UTIL_RESULTSTORAGE_H_


template <unsigned int DIM, unsigned int WIDTH>
class ResultStorage {
public:
	static const size_t CAPACITY = 1000;
	size_t nextIndex_;
	ResultStorage<DIM, WIDTH>* nextStorage_;
	Entry<DIM, WIDTH> entries_[CAPACITY];

	ResultStorage();
	~ResultStorage();
	ResultStorage<DIM, WIDTH>* add(const Entry<DIM,WIDTH>& entry);
	void recursiveDelete();
};



template <unsigned int DIM, unsigned int WIDTH>
ResultStorage<DIM, WIDTH>::ResultStorage() : nextIndex_(0), nextStorage_(NULL) {
}

template <unsigned int DIM, unsigned int WIDTH>
ResultStorage<DIM, WIDTH>::~ResultStorage() {
}

template <unsigned int DIM, unsigned int WIDTH>
ResultStorage<DIM, WIDTH>* ResultStorage<DIM, WIDTH>::add(
		const Entry<DIM, WIDTH>& entry) {
	entries_[nextIndex_++] = entry;
	if (nextIndex_ == CAPACITY) {
		nextStorage_ = new ResultStorage<DIM, WIDTH>();
		return nextStorage_;
	} else {
		return this;
	}
}


#endif /* SRC_UTIL_RESULTSTORAGE_H_ */
