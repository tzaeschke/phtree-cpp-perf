/*
 * EntryTreeMap.h
 *
 *  Created on: Aug 15, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_ENTRYTREEMAP_H_
#define SRC_UTIL_ENTRYTREEMAP_H_

template <unsigned int DIM>
class Node;
template <unsigned int DIM, unsigned int WIDTH>
class Entry;


template <unsigned int DIM, unsigned int WIDTH>
class EntryTreeMap {
public:
	EntryTreeMap();
	~EntryTreeMap();

	void clearMap();
	void enforcePreviousNode();
	void getNextUndeletedNode(Node<DIM>** lowestCommonNode, size_t* index);
	bool compareForStart(const Entry<DIM, WIDTH>& entry);
	void put(size_t startIndex, Node<DIM>* containedInNode);

	const Entry<DIM, WIDTH>* createEntry(const std::vector<unsigned long> &values, int id);

private:

	bool currentEntryIsFirst_;
	unsigned int lowestCommonNodeIndex_;
	Entry<DIM, WIDTH> entry1_;
	Entry<DIM, WIDTH> entry2_;
	size_t startIndices_[WIDTH]; // Probably a fraction of the width is also fine if using the highest value in case of exceeding
	Node<DIM>* mappedNodes_[WIDTH];

	inline const Entry<DIM, WIDTH>* currentEntry() const;
	inline const Entry<DIM, WIDTH>* lastEntry() const;
};

#include "../nodes/Node.h"
#include "../Entry.h"
#include "../util/MultiDimBitset.h"

using namespace std;

template <unsigned int DIM, unsigned int WIDTH>
EntryTreeMap<DIM, WIDTH>::EntryTreeMap() : currentEntryIsFirst_(false), lowestCommonNodeIndex_(-1u),
	entry1_(), entry2_(), startIndices_(), mappedNodes_() {
}

template <unsigned int DIM, unsigned int WIDTH>
EntryTreeMap<DIM, WIDTH>::~EntryTreeMap() {
}

template <unsigned int DIM, unsigned int WIDTH>
const Entry<DIM, WIDTH>* EntryTreeMap<DIM, WIDTH>::createEntry(const std::vector<unsigned long> &values, int id) {
	currentEntryIsFirst_ = !currentEntryIsFirst_;
	if (currentEntryIsFirst_) {
		entry1_.reinit(values, id);
		return &entry1_;
	} else {
		entry2_.reinit(values, id);
		return &entry2_;
	}
}

template <unsigned int DIM, unsigned int WIDTH>
const Entry<DIM, WIDTH>* EntryTreeMap<DIM, WIDTH>::currentEntry() const {
	if (currentEntryIsFirst_) {
		return &entry1_;
	} else {
		return &entry2_;
	}
}

template <unsigned int DIM, unsigned int WIDTH>
const Entry<DIM, WIDTH>* EntryTreeMap<DIM, WIDTH>::lastEntry() const {
	if (currentEntryIsFirst_) {
		return &entry2_;
	} else {
		return &entry1_;
	}
}

template <unsigned int DIM, unsigned int WIDTH>
void EntryTreeMap<DIM, WIDTH>::clearMap() {
	lowestCommonNodeIndex_ = -1u;
}

template <unsigned int DIM, unsigned int WIDTH>
bool EntryTreeMap<DIM, WIDTH>::compareForStart(const Entry<DIM, WIDTH>& entry) {
	if (lowestCommonNodeIndex_ != -1u) { // TODO what about index = 0?
		// compare the new entry to the old one and reset reset the old entry
		const size_t differentAtMsb = (!lastEntry())? 0
				: MultiDimBitset<DIM>::compareFullAligned(entry.values_, DIM * WIDTH, lastEntry()->values_);
		if (differentAtMsb == WIDTH) { return true; }

		// check which node the insertion needs to start at, i.e. the lowest
		// common node that does not require any changes
		unsigned int l = 0;
		unsigned int r = lowestCommonNodeIndex_;
		// find the lowest common node using binary search
		while (l < r) {
			// check interval [l, r) => never uses an index after the highest valid one
			const unsigned int middle = (l + r) / 2;
			const size_t currentWidthStart = startIndices_[middle];
			if (currentWidthStart < differentAtMsb) {
				l = middle + 1;
			} else if (currentWidthStart > differentAtMsb) {
				r = middle;
			} else {
				// exact match so need to take the previous node as the first
				// bit of the found node is different from the one needed
				l = middle;
				r = middle;
			}
		}

		// might have to go one step back
		if (startIndices_[l] > differentAtMsb) { l -= 1; }

#ifdef PRINT
		cout << "new is different at index " << differentAtMsb << endl;
		for (unsigned i = 0; i <= lowestCommonNodeIndex_; ++i) {
			cout << i << "\t-> [" << startIndices_[i] << "-";
			if (i != lowestCommonNodeIndex_) { cout << (startIndices_[i + 1] - 1) << "]"; }
			else { cout << (WIDTH - 1) << "]"; }

			if (l - 1 == i) { cout << "\tstart here"; }
			else if (l == i) { cout << "\tmatch here"; }
			cout << endl;
		}
#endif

		// the entry requires the same path from the root node to the node at
		// index l - 1 as it differs at the node at index l
		lowestCommonNodeIndex_ = (l == -1u)? -1u : l - 1;

	}

	return false;
}

template <unsigned int DIM, unsigned int WIDTH>
void EntryTreeMap<DIM, WIDTH>::getNextUndeletedNode(Node<DIM>** lowestCommonNode, size_t* index) {
	assert (lowestCommonNodeIndex_ != (-1uL));

	for (unsigned i = lowestCommonNodeIndex_; i != (-1u); --i) {
		if (!mappedNodes_[i]->removed) {
			lowestCommonNodeIndex_ = i;
			(*lowestCommonNode) = mappedNodes_[i];
			(*index) = startIndices_[i];

#ifdef PRINT
			cout << "deep start at node " << (i + 1) << ", index " << (*index) << "/" << WIDTH << " -> " << flush;
#endif

			return;
		}
	}

	lowestCommonNodeIndex_ = -1u;
	(*lowestCommonNode) = NULL;
	(*index) = 0;
}

template <unsigned int DIM, unsigned int WIDTH>
void EntryTreeMap<DIM, WIDTH>::enforcePreviousNode() {
	assert (lowestCommonNodeIndex_ != (-1uL));
	--lowestCommonNodeIndex_;
}

template <unsigned int DIM, unsigned int WIDTH>
void EntryTreeMap<DIM, WIDTH>::put(size_t startIndex, Node<DIM>* containedInNode) {
	assert (containedInNode && !containedInNode->removed);
	assert (startIndex < WIDTH);
	assert (lowestCommonNodeIndex_ == -1u || startIndices_[lowestCommonNodeIndex_] < startIndex);
	assert (lowestCommonNodeIndex_ == -1u || lowestCommonNodeIndex_ < WIDTH);
	++lowestCommonNodeIndex_;
	startIndices_[lowestCommonNodeIndex_] = startIndex;
	mappedNodes_[lowestCommonNodeIndex_] = containedInNode;
}


#endif /* SRC_UTIL_ENTRYTREEMAP_H_ */
