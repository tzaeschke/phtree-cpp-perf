/*
 * TEntryBuffer.h
 *
 *  Created on: Jul 1, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_TENTRYBUFFER_H_
#define SRC_UTIL_TENTRYBUFFER_H_

#include <assert.h>

template <unsigned int DIM>
class Node;

template <unsigned int DIM>
class TEntryBuffer {
public:
	TEntryBuffer();
	virtual ~TEntryBuffer() {};

	void updateNode(Node<DIM>* node);
	std::pair<Node<DIM>*, unsigned long> getNodeAndAddress();

protected:
	Node<DIM>* node_;
	unsigned long nodeHcAddress;
};

#include "../nodes/Node.h"

using namespace std;

template <unsigned int DIM>
TEntryBuffer<DIM>::TEntryBuffer() : node_(NULL), nodeHcAddress(0) {
}

template <unsigned int DIM>
void TEntryBuffer<DIM>::updateNode(Node<DIM>* node) {
	node_ = node;
}

template <unsigned int DIM>
pair<Node<DIM>*, unsigned long> TEntryBuffer<DIM>::getNodeAndAddress() {
	return pair<Node<DIM>*, unsigned long>(node_, nodeHcAddress);
}

#endif /* SRC_UTIL_TENTRYBUFFER_H_ */
