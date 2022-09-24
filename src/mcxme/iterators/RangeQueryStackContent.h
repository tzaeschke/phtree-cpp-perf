/*
 * RangeQueryStackContent.h
 *
 *  Created on: May 18, 2016
 *      Author: max
 */

#ifndef SRC_ITERATORS_RANGEQUERYSTACKCONTENT_H_
#define SRC_ITERATORS_RANGEQUERYSTACKCONTENT_H_

#include "../nodes/Node.h"
#include "../iterators/NodeIterator.h"

template <unsigned int DIM>
class Node;

template <unsigned int DIM>
struct RangeQueryStackContent {
public:
	bool fullyContained;
	bool lowerContained;
	bool upperContained;

	unsigned int lowerMask_;
	unsigned int upperMask_;

	unsigned int prefixLength_;

	unsigned long lowerCompEqual;
	unsigned long lowerCompSmaller;

	unsigned long upperCompEqual;
	unsigned long upperCompSmaller;

	const Node<DIM>* node_;
	// TODO make iterators lokal
	NodeIterator<DIM>* startIt_;
	NodeIterator<DIM>* endIt_;
};

#endif /* SRC_ITERATORS_RANGEQUERYSTACKCONTENT_H_ */
