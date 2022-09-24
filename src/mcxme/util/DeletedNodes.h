/*
 * DeletedNodes.h
 *
 *  Created on: Aug 15, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_DELETEDNODES_H_
#define SRC_UTIL_DELETEDNODES_H_

template <unsigned int DIM>
class Node;

template <unsigned int DIM>
class DeletedNodes {
public:
	DeletedNodes();
	~DeletedNodes();


	bool full() const;
	void deleteAll();
	void add(Node<DIM>* node);

private:

	static const size_t SIZE = 2000;

	size_t nextIndex_;
	Node<DIM>* buffer_[SIZE];
};

#include "../nodes/Node.h"

using namespace std;

template <unsigned int DIM>
DeletedNodes<DIM>::DeletedNodes() : nextIndex_(0), buffer_() {

}

template <unsigned int DIM>
DeletedNodes<DIM>::~DeletedNodes() {
	deleteAll();
}

template<unsigned int DIM>
bool DeletedNodes<DIM>::full() const {
	return nextIndex_ == SIZE;
}

template<unsigned int DIM>
void DeletedNodes<DIM>::deleteAll() {
	for (unsigned i = 0; i < nextIndex_; ++i) {
		assert (buffer_[i]);
		delete buffer_[i];
		buffer_[i] = NULL;
	}

	nextIndex_ = 0;
}

template <unsigned int DIM>
void DeletedNodes<DIM>::add(Node<DIM>* node) {
	assert (!node->removed);
	node->removed = true;
	assert (nextIndex_ < SIZE);
	assert (!buffer_[nextIndex_]);
	buffer_[nextIndex_++] = node;
}


#endif /* SRC_UTIL_DELETEDNODES_H_ */
