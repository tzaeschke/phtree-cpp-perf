/*
 * RangeQueryIterator.h
 *
 *  Created on: Mar 16, 2016
 *      Author: max
 */

#ifndef SRC_ITERATORS_RANGEQUERYITERATOR_H_
#define SRC_ITERATORS_RANGEQUERYITERATOR_H_

#include <stack>
#include "../nodes/Node.h"
#include "../iterators/NodeIterator.h"
#include "../iterators/RangeQueryStackContent.h"
#include "../util/MultiDimBitset.h"
#include "../Entry.h"

template <unsigned int DIM>
class Node;

template <unsigned int DIM, unsigned int WIDTH>
class RangeQueryIterator {
public:
	// TODO add version that does not recreate the spatial data but only returns IDs
	RangeQueryIterator(std::vector<std::pair<unsigned long, const Node<DIM>*>>* nodeStack,
			const Entry<DIM, WIDTH>& lowerLeft,
			const Entry<DIM, WIDTH>& upperRight);
	virtual ~RangeQueryIterator();

	Entry<DIM, WIDTH> next();
	bool hasNext() const;

private:
	static const size_t highestAddress = (1u << DIM) - 1;

	bool hasNext_;
	size_t currentIndex_;

	// TODO replace with array implementation
	// stack of not fully traversed nodes
	std::stack<RangeQueryStackContent<DIM>> stack_;
	// relevant information for the currently processed node
	RangeQueryStackContent<DIM> currentContent;
	// storage for bits from higher nodes (lower levels from stack)
	// Combined with a suffix this defines an entry.
	unsigned long currentValue[1 + (DIM * WIDTH - 1) / (sizeof (unsigned long) * 8)];
	// The address contents of the currently processed address in the currently processed node
	NodeAddressContent<DIM> currentAddressContent;

	const Entry<DIM, WIDTH> lowerLeftCorner_;
	const Entry<DIM, WIDTH> upperRightCorner_;


	void stepUp();
	bool stepDown(const Node<DIM>* nextNode, unsigned long hcAddress);
	inline bool isInMaskRange(unsigned long hcAddress) const;
	inline bool isSuffixInRange();
	inline void createCurrentContent(const Node<DIM>* nextNode, unsigned long hcAddress, unsigned int prefixLength);
	inline void goToNextValidSuffix();

	// utility functions for local calculations:
	inline void compareLocalLower(unsigned long prevLeft, unsigned long currLeft,
				unsigned long prevAddress,
				unsigned long* outCompEqual, unsigned long* outCompSmaller) const;
	inline void compareLocalUpper(unsigned long prevRight, unsigned long currRight,
				unsigned long prevAddress,
				unsigned long* outCompEqual, unsigned long* outCompSmaller) const;
	inline void compareLocalRangeLower(unsigned int skipLowerNBits, unsigned int compareNBits,
			unsigned long* outCompEqual, unsigned long* outCompSmaller) const;
	inline void compareLocalRangeUpper(unsigned int skipLowerNBits, unsigned int compareNBits,
			unsigned long* outCompEqual, unsigned long* outCompSmaller) const;
	inline void connectLocalToPrevLower(unsigned long prevEqualComp, unsigned long prevSmallerComp,
			unsigned long localEqualComp, unsigned long localSmallerComp,
			unsigned long* outNewEqualComp, unsigned long* outNewSmallerComp) const;
	inline void connectLocalToPrevUpper(unsigned long prevEqualComp, unsigned long prevSmallerComp,
			unsigned long localEqualComp, unsigned long localSmallerComp,
			unsigned long* outNewEqualComp, unsigned long* outNewSmallerlComp) const;
};

#include <assert.h>
#include "../nodes/NodeAddressContent.h"
#include "../util/SpatialSelectionOperationsUtil.h"

using namespace std;

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>::RangeQueryIterator(vector<pair<unsigned long, const Node<DIM>*>>* visitedNodes,
		const Entry<DIM, WIDTH>& lowerLeft, const Entry<DIM, WIDTH>& upperRight) : hasNext_(true),
		currentIndex_(0), stack_(),
		currentValue(), lowerLeftCorner_(lowerLeft),
		upperRightCorner_(upperRight) {

#ifndef NDEBUG
	// validation only: lower left < upper right
	pair<unsigned long, unsigned long> comp = MultiDimBitset<DIM>::
				compareSmallerEqual(lowerLeftCorner_.values_, upperRightCorner_.values_, DIM * WIDTH, 0, highestAddress);
	assert ((comp.first == highestAddress) && "should be: lower left < upper right");
#endif

	assert (!visitedNodes->empty() && "at least the root node must have been visited");
	if (visitedNodes->empty()) {
		hasNext_ = false;
	} else {
		currentContent.fullyContained = false;
		currentContent.lowerContained = false;
		currentContent.upperContained = false;
		// start: range window is fully included in domain
		currentContent.lowerCompEqual = highestAddress;
		currentContent.lowerCompSmaller = 0;
		currentContent.upperCompEqual = highestAddress;
		currentContent.upperCompSmaller = 0;
		// the first node has to be the root node which does not have a prefix!
		const Node<DIM>* root = (*visitedNodes)[0].second;
		createCurrentContent(root, 0, 0);
		for (unsigned int i = 1; i < visitedNodes->size(); ++i) {
			const pair<unsigned long, const Node<DIM>*> nextNode = (*visitedNodes)[i];
			// TODO actually no need to validate prefixes since the lookup already did that?!
			bool prefixIncluded = stepDown(nextNode.second, nextNode.first);
			assert (prefixIncluded);
		}

		// thereby checks if there is any valid entry in the range
		goToNextValidSuffix();
	}
}

template <unsigned int DIM, unsigned int WIDTH>
RangeQueryIterator<DIM, WIDTH>::~RangeQueryIterator() {
	while (!stack_.empty()) {
		RangeQueryStackContent<DIM> content = stack_.top();
		delete content.startIt_;
		delete content.endIt_;
		stack_.pop();
	}

	if (currentContent.startIt_)
		delete currentContent.startIt_;
	if (currentContent.endIt_)
		delete currentContent.endIt_;
}

template <unsigned int DIM, unsigned int WIDTH>
Entry<DIM, WIDTH> RangeQueryIterator<DIM, WIDTH>::next() {
	assert (hasNext());
	assert (currentAddressContent.exists && isInMaskRange(currentAddressContent.address)
		&& !currentAddressContent.hasSubnode && isSuffixInRange());
	assert (currentIndex_ < WIDTH);
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM, 0)
			&& "there always needs to be space for the last interleaved address");
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));

	// found a valid suffix in the range
	Entry<DIM, WIDTH> entry(currentValue, currentAddressContent.id);
	// copy the suffix into the entry
	MultiDimBitset<DIM>::pushBackValue(currentAddressContent.address, entry.values_, DIM * (WIDTH - currentIndex_ - 1));
	if (currentIndex_ < WIDTH - 1) {
		// need to copy the remaining suffix into the block
		const size_t suffixBits = DIM * (WIDTH - (currentIndex_ + 1));
		assert (suffixBits > 0);
		assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, suffixBits, 0));
		MultiDimBitset<DIM>::pushBackBitset(currentAddressContent.getSuffixStartBlock(), suffixBits, entry.values_, 0);
	}


#ifndef NDEBUG
	// validation only: is the entry contained in the current node
	assert (entry.id_ == currentContent.node_->lookup(currentAddressContent.address, true).id);

	// validation only: is the retrieved entry part of the tree?
	const Node<DIM>* rootNode = NULL;
	if (stack_.empty()) rootNode = currentContent.node_;
	else {
		// the root node is at the bottom of the stack
		stack<RangeQueryStackContent<DIM>> reverseStack;
		while (!stack_.empty()) {
			RangeQueryStackContent<DIM> top = stack_.top();
			rootNode = top.node_;
			reverseStack.push(top);
			stack_.pop();
		}
		while (!reverseStack.empty()) {
			stack_.push(reverseStack.top());
			reverseStack.pop();
		}
	}

	std::pair<bool, int> lookup = SpatialSelectionOperationsUtil<DIM, WIDTH>::lookup(entry, rootNode, NULL);
	assert (lookup.first && lookup.second == entry.id_);

	// validation only: lower left <= entry <= upper right
	// i.e. if in any dimension the inverse operation for < is not 0 there is an error
	pair<unsigned long, unsigned long> lowerComp = MultiDimBitset<DIM>::
				compareSmallerEqual(lowerLeftCorner_.values_, entry.values_, DIM * WIDTH, 0, highestAddress);
	pair<unsigned long, unsigned long> upperComp = MultiDimBitset<DIM>::
				compareSmallerEqual(entry.values_, upperRightCorner_.values_, DIM * WIDTH, 0, highestAddress);
	assert (((lowerComp.first | lowerComp.second) == highestAddress) && "should be: lower left <= entry");
	assert (((upperComp.first | upperComp.second) == highestAddress) && "should be: entry <= upper right");
#endif

	++(*currentContent.startIt_);
	goToNextValidSuffix();

	return entry;
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::goToNextValidSuffix() {
	assert (hasNext_);

	do {
		while ((*currentContent.startIt_) == (*currentContent.endIt_) && hasNext_) {
			// ascend to a previous level if the end of the node was reached
			stepUp();
			++(*currentContent.startIt_);
			assert (!hasNext_ || (*currentContent.startIt_) <= (*currentContent.endIt_));
		}

		if (!hasNext_) break;
		currentAddressContent = *(*currentContent.startIt_);
		assert (currentAddressContent.exists);
		if (!isInMaskRange(currentAddressContent.address)) {
			++(*currentContent.startIt_);
		} else if (currentAddressContent.hasSubnode) {
			// descend to the next level in case of a subnode
			bool prefixIncluded = stepDown(currentAddressContent.subnode, currentAddressContent.address);
			if (!prefixIncluded) {
				// the prefix of the current node was not included so do not descend
				++(*currentContent.startIt_);
			}
		} else if (isSuffixInRange()) {
			// found a suffix with a valid address
			break;
		} else {
			// the suffix was invalid so continue the search
			++(*currentContent.startIt_);
		}
	} while (hasNext_);

	assert (!hasNext_ || (currentAddressContent.exists && isInMaskRange(currentAddressContent.address) && !currentAddressContent.hasSubnode));
}

template <unsigned int DIM, unsigned int WIDTH>
bool RangeQueryIterator<DIM, WIDTH>::hasNext() const {
	return hasNext_;
}

template <unsigned int DIM, unsigned int WIDTH>
bool RangeQueryIterator<DIM, WIDTH>::isInMaskRange(unsigned long hcAddress) const {

	assert (currentContent.upperMask_ < (1uL << DIM));
	assert (currentContent.lowerMask_ <= currentContent.upperMask_);
	assert (currentContent.lowerMask_ <= hcAddress && hcAddress <= currentContent.upperMask_);

	const bool addressMatch = currentContent.fullyContained
			|| (((hcAddress | currentContent.lowerMask_) & currentContent.upperMask_) == hcAddress);
	return addressMatch;
}

template <unsigned int DIM, unsigned int WIDTH>
bool RangeQueryIterator<DIM, WIDTH>::isSuffixInRange() {
	assert (currentAddressContent.exists && !currentAddressContent.hasSubnode
			&& isInMaskRange(currentAddressContent.address));
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));

	if (currentContent.fullyContained) { return true; }

	// <prev. compared>
	//                 <---- local comparison --->
	// [    higher    |curr|       suffix        ]

	// Verify if the final entry drops out of the range!
	// assemble the entry from the buffer, the current HC address and the current suffix
	const unsigned int suffixLength = WIDTH - currentIndex_ - 1;
	const unsigned int suffixBits = DIM * suffixLength;
	MultiDimBitset<DIM>::pushBackValue(currentAddressContent.address, currentValue, DIM * (WIDTH - currentIndex_ - 1));
	if (suffixLength > 0)
		MultiDimBitset<DIM>::pushBackBitset(currentAddressContent.getSuffixStartBlock(), suffixBits, currentValue, 0);

	bool suffixContained = true;
	// local comparison: compare the entire suffix and the current address
	const unsigned int compareLowestNBits = suffixBits + DIM;
	// validate: range lower left <= entry (suffix)
	if (!currentContent.lowerContained) {
		unsigned long lowerCompRangeEqualLocal;
		unsigned long lowerCompRangeSmallerLocal;
		compareLocalRangeLower(0, compareLowestNBits, &lowerCompRangeEqualLocal, &lowerCompRangeSmallerLocal);
		unsigned long lowerCompEqualWithSuffix;
		unsigned long lowerCompSmallerWithSuffix;
		connectLocalToPrevLower(currentContent.lowerCompEqual, currentContent.lowerCompSmaller,
				lowerCompRangeEqualLocal, lowerCompRangeSmallerLocal,
				&lowerCompEqualWithSuffix, &lowerCompSmallerWithSuffix);

		// lower <= current (with prefix)
		suffixContained = (lowerCompEqualWithSuffix | lowerCompSmallerWithSuffix) == highestAddress;
	}

	// validate: entry (suffix) <= range upper right
	if (suffixContained && !currentContent.upperContained) {
		unsigned long upperCompRangeEqualLocal;
		unsigned long upperCompRangeSmallerLocal;
		compareLocalRangeUpper(0, compareLowestNBits, &upperCompRangeEqualLocal, &upperCompRangeSmallerLocal);
		unsigned long upperCompEqualWithSuffix;
		unsigned long upperCompSmallerWithSuffix;
		connectLocalToPrevUpper(currentContent.upperCompEqual, currentContent.upperCompSmaller,
				upperCompRangeEqualLocal, upperCompRangeSmallerLocal,
				&upperCompEqualWithSuffix, &upperCompSmallerWithSuffix);

		// upper >= current (with prefix) <=> (not(smaller) and not(equal)) or equal <=> not(smaller) or equal
		suffixContained = (upperCompEqualWithSuffix | (highestAddress & (~upperCompSmallerWithSuffix))) == highestAddress;
	}

	// clear the buffer (remove the HC address and possibly also the suffix)
	if (suffixLength > 0) {
		MultiDimBitset<DIM>::removeHighestBits(currentValue, (1 + suffixLength) * DIM, (1 + suffixLength) * DIM);
	} else {
		MultiDimBitset<DIM>::clearValue(currentValue, DIM * (WIDTH - currentIndex_ - 1));
	}

	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
	if (suffixContained) {
		// the suffix is in the upper range
		// thus it will be in the next entry to be returned
		// TODO --> no need to clear the value
		return true;
	} else {
		return false;
	}
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::stepUp() {
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
	assert ((*currentContent.startIt_) == (*currentContent.endIt_));
	// TODO if several stack contents are skipped only a single remove operation is needed!
	if (stack_.empty()) {
		hasNext_ = false;
		assert (currentContent.prefixLength_ == 0
				&& "the last node should be the root which does not have a prefix");
		assert (currentIndex_ == 0);
	} else {
		// remove the prefix of the last node
		const size_t prefixLength = currentContent.prefixLength_;
		if (prefixLength > 0) {
			// remove current interleaved address (1) and prefix (+ prefixLength) bits (*DIM)
			currentIndex_ -= (prefixLength + 1);
			const unsigned int freeLsbBits = DIM * (WIDTH - currentIndex_);
			MultiDimBitset<DIM>::removeHighestBits(currentValue, freeLsbBits, (prefixLength + 1) * DIM);
			assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
		} else {
			// only need to clear the current interleaved address
			MultiDimBitset<DIM>::clearValue(currentValue, DIM * (WIDTH - currentIndex_));
			currentIndex_ -= 1;
			assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
		}

		// clear memory
		delete currentContent.startIt_;
		delete currentContent.endIt_;

		// restore the last contents from the top of the stack
		currentContent = stack_.top();
		stack_.pop();
	}
}

template <unsigned int DIM, unsigned int WIDTH>
bool RangeQueryIterator<DIM, WIDTH>::stepDown(const Node<DIM>* nextNode, unsigned long hcAddress) {
	assert (nextNode && hcAddress < (1uL << DIM));
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
	assert ((*currentContent.startIt_) < (*currentContent.endIt_));

	// msb               [interleaved format]                      lsb
	// <-filled-><DIM><DIM *|prefix|><DIM><-------- ignored --------->
	// <prev. compared>             <------- later compared --------->
	//          < local comparison >
	// [ higher |last|current prefix|curr|     lower node bits       ]
	// last: node currently descending from
	// curr: node currently descending into

	// add the current interleaved address
	currentIndex_ += 1;
	MultiDimBitset<DIM>::pushBackValue(hcAddress, currentValue, (WIDTH - currentIndex_) * DIM);
	const size_t prefixLength = nextNode->getPrefixLength();
	unsigned long lowerCompEqualWithPrefix = currentContent.lowerCompEqual;
	unsigned long lowerCompSmallerWithPrefix = currentContent.lowerCompSmaller;
	unsigned long upperCompEqualWithPrefix = currentContent.upperCompEqual;
	unsigned long upperCompSmallerWithPrefix = currentContent.upperCompSmaller;
	unsigned long lastAddress = hcAddress;
	if (prefixLength > 0) {
		// add the prefix of the next node to the current prefix
		currentIndex_ += prefixLength;
		const unsigned int ignoreNLowestBits = DIM * (WIDTH - currentIndex_);
		MultiDimBitset<DIM>::pushBackBitset(nextNode->getFixPrefixStartBlock(), prefixLength * DIM,
				currentValue, (WIDTH - currentIndex_) * DIM);
		assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));

		const unsigned int compareNBits = DIM * prefixLength + DIM;
		// varify if the prefix is still within the range
		bool prefixValid = true;
		if (!currentContent.lowerContained) {
			unsigned long lowerCompRangeEqualLocal;
			unsigned long lowerCompRangeSmallerLocal;
			compareLocalRangeLower(ignoreNLowestBits, compareNBits, &lowerCompRangeEqualLocal, &lowerCompRangeSmallerLocal);
			connectLocalToPrevLower(currentContent.lowerCompEqual, currentContent.lowerCompSmaller,
					lowerCompRangeEqualLocal, lowerCompRangeSmallerLocal,
					&lowerCompEqualWithPrefix, &lowerCompSmallerWithPrefix);

			// lower <= current (with prefix)
			prefixValid = (lowerCompEqualWithPrefix | lowerCompSmallerWithPrefix) == highestAddress;
		}

		if (prefixValid && !currentContent.upperContained) {
			unsigned long upperCompRangeEqualLocal;
			unsigned long upperCompRangeSmallerLocal;
			compareLocalRangeUpper(ignoreNLowestBits, compareNBits, &upperCompRangeEqualLocal, &upperCompRangeSmallerLocal);
			connectLocalToPrevUpper(currentContent.upperCompEqual, currentContent.upperCompSmaller,
					upperCompRangeEqualLocal, upperCompRangeSmallerLocal,
					&upperCompEqualWithPrefix, &upperCompSmallerWithPrefix);

			// upper >= current (with prefix) <=> (not(smaller) and not(equal)) or equal <=> not(smaller) or equal
			prefixValid = (upperCompEqualWithPrefix | (highestAddress & (~upperCompSmallerWithPrefix))) == highestAddress;
		}

		if (!prefixValid) {
			// the prefix was not in the range so remove it
			currentIndex_ -= (prefixLength + 1);
			const unsigned int freeLsbBits = DIM * (WIDTH - currentIndex_);
			MultiDimBitset<DIM>::removeHighestBits(currentValue, freeLsbBits, (prefixLength + 1) * DIM);
			assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
			return false;
		} else {
			lastAddress = MultiDimBitset<DIM>::interleaveBits(currentValue, currentIndex_ - 1, DIM * WIDTH);
		}

	#ifndef NDEBUG
		// validation: is the locally composed comparison for the prefix correct so far?
		pair<unsigned long, unsigned long> fullLowerComp = MultiDimBitset<DIM>::
				compareSmallerEqual(lowerLeftCorner_.values_, currentValue, DIM * WIDTH, ignoreNLowestBits, highestAddress);
		assert (fullLowerComp.first == lowerCompSmallerWithPrefix && fullLowerComp.second == lowerCompEqualWithPrefix);
		pair<unsigned long, unsigned long> fullUpperComp = MultiDimBitset<DIM>::
				compareSmallerEqual(upperRightCorner_.values_, currentValue, DIM * WIDTH, ignoreNLowestBits, highestAddress);
		assert (fullUpperComp.first == upperCompSmallerWithPrefix && fullUpperComp.second == upperCompEqualWithPrefix);
	#endif
	}

	// puts a duplicate on the stack
	stack_.push(currentContent);
	currentContent.lowerCompEqual = lowerCompEqualWithPrefix;
	currentContent.lowerCompSmaller = lowerCompSmallerWithPrefix;
	currentContent.upperCompEqual = upperCompEqualWithPrefix;
	currentContent.upperCompSmaller = upperCompSmallerWithPrefix;

	createCurrentContent(nextNode, lastAddress, prefixLength);
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, DIM * (WIDTH - currentIndex_), 0));
	return true;
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::compareLocalLower(unsigned long prevLeft, unsigned long currLeft,
			unsigned long prevAddress,
			unsigned long* outCompEqual, unsigned long* outCompSmaller) const {
	// local lower equal if: (previous lower == previous address) and (current lower == 0)
	(*outCompEqual) = (~(prevLeft ^ prevAddress)) & (~currLeft) & highestAddress;
	// local lower smaller if: previous lower (0) <  previous address (1)
	(*outCompSmaller) = (~prevLeft) & prevAddress;
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::compareLocalUpper(unsigned long prevRight, unsigned long currRight,
			unsigned long prevAddress,
			unsigned long* outCompEqual, unsigned long* outCompSmaller) const {
	// local upper | equal if: (previous upper == previous address) and (current upper == 1)
	(*outCompEqual) = (~(prevRight ^ prevAddress)) & currRight;
	// local upper | smaller in cases (u - upper range, x - current node):
	// [u|_|x|x] or [_|u|x|x] or [_|_|xu|x] or [xu|x|_|_]
	(*outCompSmaller) = highestAddress &
			(((~prevRight) & prevAddress)
			| (prevRight & prevAddress & (~currRight))
			| ((~prevRight) & (~prevAddress) & (~currRight)));
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::connectLocalToPrevLower(unsigned long prevEqualComp, unsigned long prevSmallerComp,
		unsigned long localEqualComp, unsigned long localSmallerComp,
		unsigned long* outNewEqualComp, unsigned long* outNewSmallerComp) const {
	// it is equal if it was equal and it is locally equal or if it was bigger and is locally equal
	(*outNewEqualComp) = localEqualComp & (prevEqualComp | (~prevSmallerComp));
	// it is smaller if it was smaller or if it is locally smaller
	(*outNewSmallerComp) = localSmallerComp | prevSmallerComp;
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::connectLocalToPrevUpper(unsigned long prevEqualComp, unsigned long prevSmallerComp,
		unsigned long localEqualComp, unsigned long localSmallerComp,
		unsigned long* outNewEqualComp, unsigned long* outNewSmallerComp) const {
	// it is equal if it was equal or smaller and is now equal
	(*outNewEqualComp) = localEqualComp & (prevEqualComp | prevSmallerComp);
	// it is smaller if it was equal or smaller and is now smaller
	(*outNewSmallerComp) = localSmallerComp & (prevEqualComp | prevSmallerComp);
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::compareLocalRangeLower(unsigned int skipLowerNBits, unsigned int compareNBits,
			unsigned long* outCompEqual, unsigned long* outCompSmaller) const {
	// all dimensions need to be compared where the lower left range is not already < current
	unsigned long compareDimMask = currentContent.lowerCompEqual | (highestAddress & (~currentContent.lowerCompSmaller));
	pair<unsigned long, unsigned long> lowerComp = MultiDimBitset<DIM>::compareSmallerEqual(
			lowerLeftCorner_.values_, currentValue,
			skipLowerNBits + compareNBits, skipLowerNBits,
			compareDimMask);
	(*outCompEqual) = lowerComp.second;
	(*outCompSmaller) = lowerComp.first;
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::compareLocalRangeUpper(unsigned int skipLowerNBits, unsigned int compareNBits,
			unsigned long* outCompEqual, unsigned long* outCompSmaller) const {
	// all dimensions need to be compared where the upper right range is not already > current
	unsigned long compareDimMask = currentContent.upperCompEqual | currentContent.upperCompSmaller;
	pair<unsigned long, unsigned long> upperComp = MultiDimBitset<DIM>::compareSmallerEqual(
			upperRightCorner_.values_, currentValue,
			skipLowerNBits + compareNBits, skipLowerNBits,
			compareDimMask);
	(*outCompEqual) = upperComp.second;
	(*outCompSmaller) = upperComp.first;
}

template <unsigned int DIM, unsigned int WIDTH>
void RangeQueryIterator<DIM, WIDTH>::createCurrentContent(const Node<DIM>* nextNode,
		unsigned long hcAddress, unsigned int prefixLength) {

	currentContent.node_ = nextNode;
	currentContent.prefixLength_ = prefixLength;
	assert (nextNode->getPrefixLength() == prefixLength);

#ifndef NDEBUG
	const unsigned int ignoreNLowestBits = DIM * (WIDTH - currentIndex_ - 1);
	// assert (ignoreNLowestBits >= 0 && ignoreNLowestBits <= WIDTH * DIM);
    assert (ignoreNLowestBits <= WIDTH * DIM); // TODO TZ warning fixed
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, ignoreNLowestBits + DIM, 0));
	// validate that the node (with the prefix) is within the range
	if (currentIndex_ > 0) {
		pair<unsigned long, unsigned long> checkLowerComp = MultiDimBitset<DIM>::
					compareSmallerEqual(lowerLeftCorner_.values_, currentValue, DIM * WIDTH, ignoreNLowestBits + DIM, highestAddress);
		pair<unsigned long, unsigned long> checkUpperComp = MultiDimBitset<DIM>::
					compareSmallerEqual(currentValue, upperRightCorner_.values_, DIM * WIDTH, ignoreNLowestBits + DIM, highestAddress);
		assert (((checkLowerComp.first | checkLowerComp.second) == highestAddress) && "should be: lower left <= entry");
		assert (((checkUpperComp.first | checkUpperComp.second) == highestAddress) && "should be: entry <= upper right");
	}
#endif

	// msb       [interleaved format]          lsb
	// <-filled-><DIM><DIM *|prefix|><DIM><-------- ignored --------->
	// <----- already compared -----><loc>
	// [ higher |last|current prefix|curr|     lower node bits       ]
	// last: node currently descending from
	// curr: node currently descending into
	// loc:  do a local comparison of the current value and the previous
	//       value (might be in prefix) with the proper range values

	// calculate current local lower comparisons
	unsigned long lowerCompEqualLocal = highestAddress;
	unsigned long lowerCompSmallerLocal = 0;
	if (!currentContent.fullyContained && !currentContent.lowerContained) {
		// <-filled-><DIM><-------- ignored --------->
		// [ higher |00000|     lower node bits      ]
		const unsigned long prevLowerHC = (currentIndex_ == 0)? 0 : MultiDimBitset<DIM>::interleaveBits(lowerLeftCorner_.values_, currentIndex_ - 1, DIM * WIDTH);
		const unsigned long lowerHC = MultiDimBitset<DIM>::interleaveBits(lowerLeftCorner_.values_, currentIndex_, DIM * WIDTH);
		compareLocalLower(prevLowerHC, lowerHC, hcAddress, &lowerCompEqualLocal, &lowerCompSmallerLocal);
	}

	// combine lower local comparison with previous lower comparison
	if (!currentContent.lowerContained) {
		const unsigned long lastLowerCompEqual = currentContent.lowerCompEqual;
		const unsigned long lastLowerCompSmaller = currentContent.lowerCompSmaller;
		connectLocalToPrevLower(lastLowerCompEqual, lastLowerCompSmaller, lowerCompEqualLocal, lowerCompSmallerLocal,
				&(currentContent.lowerCompEqual), &(currentContent.lowerCompSmaller));

		// lower mask: for i=0 to DIM - 1 do
		// 		lowerMask[i] = 0 <=> current value (in dimension i) is less or equal to the lower left corner (in dimension i)
		currentContent.lowerMask_ =  highestAddress & (~(currentContent.lowerCompSmaller | currentContent.lowerCompEqual));
		currentContent.lowerContained = (currentContent.lowerCompEqual == 0)
				&& (currentContent.lowerCompSmaller == highestAddress);
	}

	// set the start iterator
	if (currentContent.lowerContained) { // TODO or if the value is 0
		currentContent.startIt_ = nextNode->begin();
	} else {
		currentContent.startIt_ = nextNode->it(currentContent.lowerMask_);
	}

	// calculate current local upper comparisons
	unsigned long upperCompEqualLocal = highestAddress;
	unsigned long upperCompSmallerLocal = 0;
	if (!currentContent.fullyContained && !currentContent.upperContained) {
		// <-filled-><DIM><-------- ignored --------->
		// [ higher |11111|     lower node bits      ]
		const unsigned long prevUpperHC = (currentIndex_ == 0)? 0 : MultiDimBitset<DIM>::interleaveBits(upperRightCorner_.values_, currentIndex_ - 1, DIM * WIDTH);
		const unsigned long upperHC = MultiDimBitset<DIM>::interleaveBits(upperRightCorner_.values_, currentIndex_, DIM * WIDTH);
		compareLocalUpper(prevUpperHC, upperHC, hcAddress, &upperCompEqualLocal, &upperCompSmallerLocal);
	}

	// combine lower local comparison with previous lower comparison
	if (!currentContent.upperContained) {
		const unsigned long lastUpperCompEqual = currentContent.upperCompEqual;
		const unsigned long lastUpperCompSmaller = currentContent.upperCompSmaller;
		connectLocalToPrevUpper(lastUpperCompEqual, lastUpperCompSmaller, upperCompEqualLocal, upperCompSmallerLocal,
				&(currentContent.upperCompEqual), &(currentContent.upperCompSmaller));

		// upper mask: for i=0 to DIM - 1 do
		// 		lowerMask[i] = 1 <=> current value (in dimension i) is higher or equal to upper right corner (in dimension i)
		currentContent.upperMask_ = highestAddress & ((~currentContent.upperCompSmaller) | currentContent.upperCompEqual);
		currentContent.upperContained = (currentContent.upperCompEqual == 0)
				&& (currentContent.upperCompSmaller == 0);
	}

	// set the start iterator
	if (currentContent.upperContained) { // TODO or if the value is highestAddress
		currentContent.endIt_ = nextNode->end();
	} else {
		currentContent.endIt_ = nextNode->it(currentContent.upperMask_ + 1);
	}

	currentContent.fullyContained = currentContent.lowerMask_ == 0
					&& currentContent.upperMask_ == highestAddress
					&& currentContent.upperCompEqual == 0
					&& currentContent.lowerCompEqual == 0;

#ifndef NDEBUG
	// validation:
	// - are the local and connected values possible?
	assert ((lowerCompEqualLocal & lowerCompSmallerLocal) == 0);
	assert ((upperCompEqualLocal & upperCompSmallerLocal) == 0);
	assert ((currentContent.lowerCompEqual & currentContent.lowerCompSmaller) == 0);
	assert ((currentContent.upperCompEqual & currentContent.upperCompSmaller) == 0);

	// - are the <= comparisons for the lower range correctly assembled?
	pair<unsigned long, unsigned long> fullLowerComp = MultiDimBitset<DIM>::
				compareSmallerEqual(lowerLeftCorner_.values_, currentValue, DIM * WIDTH, ignoreNLowestBits, highestAddress);
	unsigned int fullLowerMask =  highestAddress & (~(fullLowerComp.first | fullLowerComp.second));
	assert (fullLowerComp.first == currentContent.lowerCompSmaller && fullLowerComp.second == currentContent.lowerCompEqual);
	assert (fullLowerMask == currentContent.lowerMask_);

	// - are the <= comparisons for the upper range correctly assembled?
	MultiDimBitset<DIM>::pushBackValue(highestAddress, currentValue, ignoreNLowestBits);
	pair<unsigned long, unsigned long> fullUpperComp = MultiDimBitset<DIM>::
			compareSmallerEqual(upperRightCorner_.values_, currentValue, DIM * WIDTH, ignoreNLowestBits, highestAddress);
	MultiDimBitset<DIM>::clearValue(currentValue, ignoreNLowestBits);
	assert (MultiDimBitset<DIM>::checkRangeUnset(currentValue, ignoreNLowestBits, 0));
	unsigned int fullUpperMask = highestAddress & ((~fullUpperComp.first) | fullUpperComp.second);
	assert (fullUpperComp.first == currentContent.upperCompSmaller && fullUpperComp.second == currentContent.upperCompEqual);
	assert (fullUpperMask == currentContent.upperMask_);

	// - do the mask values have the correct order and are the iterators correctly set?
	assert (currentContent.lowerMask_ <= currentContent.upperMask_ && currentContent.upperMask_ < (1uL << DIM));
	assert ((*currentContent.startIt_) <= (*currentContent.endIt_));
	assert (currentContent.lowerMask_ <= currentContent.startIt_->getAddress());
	assert (currentContent.upperMask_ < currentContent.endIt_->getAddress());

	// - if the lower and upper range is contained than the node definitely fully contained
	// notice that the inversion is not always tree
	assert (!(currentContent.lowerContained && currentContent.upperContained) || currentContent.fullyContained);
#endif
}

#endif /* SRC_ITERATORS_RANGEQUERYITERATOR_H_ */
