/*
 * EntryBuffer.h
 *
 *  Created on: Jun 30, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_ENTRYBUFFER_H_
#define SRC_UTIL_ENTRYBUFFER_H_

#include <set>
#include "../Entry.h"
#include "../util/TEntryBuffer.h"
#include <atomic>

template <unsigned int DIM>
class Node;
template <unsigned int DIM, unsigned int WIDTH>
class PHTree;
template <unsigned int DIM, unsigned int WIDTH>
class EntryBufferPool;

template <unsigned int DIM, unsigned int WIDTH>
class EntryBuffer : public TEntryBuffer<DIM> {
	friend class EntryBufferPool<DIM, WIDTH>;
public:
	EntryBuffer();
	~EntryBuffer() {};

	bool insert(const Entry<DIM, WIDTH>& entry);
	bool full() const;
	bool empty() const;
	size_t capacity() const;
	void clear();
	Entry<DIM, WIDTH>* init(size_t suffixLength, Node<DIM>* node, unsigned long hcAddress);
	void updateNode(Node<DIM>* node);
	std::pair<Node<DIM>*, unsigned long> getNodeAndAddress();
	Node<DIM>* flushToSubtree(); // TODO should be made const
	EntryBufferPool<DIM, WIDTH>* getPool();

	bool assertCleared() const;

private:

	static const size_t capacity_ = 10;

	atomic<bool> flushing_;
	atomic<size_t> nextIndex_;
	bool inUse;
	size_t suffixBits_;
	EntryBufferPool<DIM, WIDTH>* pool_;
	// - symmetric matrix: need to store n/2 (n+1) fields only
	// - stores the diagonal too for easier access
	// - stores the lower matrix because LCP values of one row are stored together
	unsigned int lcps_[capacity_ * (capacity_ + 1) / 2];
	Entry<DIM, WIDTH> buffer_[capacity_];
	bool insertCompleted_[capacity_];

	// TODO validation only:
	const Entry<DIM, WIDTH>* originals_[capacity_];

	inline void setLcp(unsigned int row, unsigned int column, unsigned int lcp);
	inline unsigned int getLcp(unsigned int row, unsigned int column) const;
	void setPool(EntryBufferPool<DIM, WIDTH>* pool);
};

#include <assert.h>
#include "../util/MultiDimBitset.h"
#include "../util/NodeTypeUtil.h"
#include "../nodes/Node.h"
#include "../PHTree.h"
#include "../util/EntryBufferPool.h"

template <unsigned int DIM, unsigned int WIDTH>
EntryBuffer<DIM, WIDTH>::EntryBuffer() : flushing_(false), nextIndex_(0), inUse(false), suffixBits_(0), pool_(NULL), lcps_(), buffer_(), insertCompleted_() {
}

template <unsigned int DIM, unsigned int WIDTH>
Entry<DIM, WIDTH>* EntryBuffer<DIM, WIDTH>::init(size_t suffixLength, Node<DIM>* node, unsigned long hcAddress) {
	assert (0 < suffixLength && suffixLength <= (WIDTH - 1));
	assert (!flushing_);
	suffixBits_ = suffixLength * DIM;
	this->node_ = node;
	this->nodeHcAddress = hcAddress;
	nextIndex_ = 1;
	flushing_ = false;
	assert (!insertCompleted_[0]);
	insertCompleted_[0] = true;
	return &(buffer_[0]);
}

template <unsigned int DIM, unsigned int WIDTH>
void EntryBuffer<DIM, WIDTH>::clear() {
	assert (suffixBits_ > 0);
	assert (flushing_);

	const size_t i = nextIndex_;
	const size_t c = capacity_;
	const size_t n = c; // TODO less work possible: min(i, c);
	const size_t blocksPerEntry = 1 + (suffixBits_ - 1) / MultiDimBitset<DIM>::bitsPerBlock;
	for (unsigned row = 0; row < n; ++row) {
		insertCompleted_[row] = false;

		for (unsigned column = 0; column <= row; ++column) {
			setLcp(row, column, 0);
		}

		for (unsigned block = 0; block < blocksPerEntry; ++block) {
			buffer_[row].values_[block] = 0;
		}
	}

	suffixBits_ = 0;
	nextIndex_ = 0;
	this->node_ = NULL;
	flushing_ = false;
}

template <unsigned int DIM, unsigned int WIDTH>
bool EntryBuffer<DIM, WIDTH>::assertCleared() const {
	assert (suffixBits_ == 0);
	assert (nextIndex_ == 0);
	assert (!flushing_);
	for (unsigned i = 0; i < capacity_; ++i) {
		assert (!insertCompleted_[i]);

		for (unsigned j = 0; j < capacity_; j++) {
			assert (getLcp(i, j) == 0);
		}
	}

	assert (!this->node_);

	return true;
}

template <unsigned int DIM, unsigned int WIDTH>
inline void EntryBuffer<DIM, WIDTH>::setLcp(unsigned int row, unsigned int column, unsigned int lcp) {
	assert (row < capacity_ && column < capacity_);
	// only stores lower half of the symmetrical matrix
	const unsigned int i = (row > column)? row : column;
	const unsigned int j = (row > column)? column : row;
	const unsigned int index = i * (i + 1) / 2 + j;
	lcps_[index] = lcp;
}

template <unsigned int DIM, unsigned int WIDTH>
inline unsigned int EntryBuffer<DIM, WIDTH>::getLcp(unsigned int row, unsigned int column) const {
	assert (row < capacity_ && column < capacity_);
	// only stores lower half of the symmetrical matrix
	const unsigned int i = (row > column)? row : column;
	const unsigned int j = (row > column)? column : row;
	const unsigned int index = i * (i + 1) / 2 + j;
	return lcps_[index];
}

template <unsigned int DIM, unsigned int WIDTH>
EntryBufferPool<DIM,WIDTH>* EntryBuffer<DIM, WIDTH>::getPool() {
	return pool_;
}

template <unsigned int DIM, unsigned int WIDTH>
void EntryBuffer<DIM, WIDTH>::setPool(EntryBufferPool<DIM,WIDTH>* pool) {
	pool_ = pool;
}

template <unsigned int DIM, unsigned int WIDTH>
size_t EntryBuffer<DIM, WIDTH>::capacity() const {
	return capacity_;
}

template <unsigned int DIM, unsigned int WIDTH>
bool EntryBuffer<DIM, WIDTH>::full() const {
	return nextIndex_ >= capacity_  && suffixBits_ > 0; // TODO replace with inUse variable from recent version;
}

template <unsigned int DIM, unsigned int WIDTH>
bool EntryBuffer<DIM, WIDTH>::empty() const {
	return nextIndex_ == 0;
}

template <unsigned int DIM, unsigned int WIDTH>
bool EntryBuffer<DIM, WIDTH>::insert(const Entry<DIM, WIDTH>& entry) {
	assert (suffixBits_ > 0 && suffixBits_ % DIM == 0);
	assert (!flushing_);
	const size_t i = nextIndex_++;
	if (i >= capacity_) { return false; }
	// exclusively responsible for index i

	assert (0 < i && i < capacity_);
	assert (!insertCompleted_[i]);
	// copy ID and necessary bits into the local buffer
	MultiDimBitset<DIM>::duplicateLowestBitsAligned(entry.values_, suffixBits_, buffer_[i].values_);
	insertCompleted_[i] = true; // TODO place after duplication?
	buffer_[i].id_ = entry.id_;
	originals_[i] = &entry;

	// compare the new entry to all previously inserted entries
	const unsigned int startIndexDim = WIDTH - (suffixBits_ / DIM);
	for (unsigned other = 0; other < i; ++other) {
		while (!insertCompleted_[other]) {}; // spin until the previous thread is done
		// TODO no need to compare all values!
		// TODO no need to compare to full values!
		assert (getLcp(other, i) == 0 && getLcp(i, other) == 0);
		assert (MultiDimBitset<DIM>::checkRangeUnset(buffer_[other].values_, DIM * WIDTH, suffixBits_));
		const pair<bool, size_t> comp = MultiDimBitset<DIM>::compare(
				entry.values_, DIM * WIDTH, startIndexDim, WIDTH, buffer_[other].values_, suffixBits_);
		assert (!flushing_);
		setLcp(i, other, comp.second);
		if (comp.first) { setLcp(i, i, -1u); }
	}

	return true;
}


template <unsigned int DIM, unsigned int WIDTH>
Node<DIM>* EntryBuffer<DIM, WIDTH>::flushToSubtree() {
	assert (suffixBits_ > 0);
	assert ((this->node_->lookup(this->nodeHcAddress, true).exists)
				&& (this->node_->lookup(this->nodeHcAddress, true).hasSpecialPointer));

	// builds the subtree from bottom up
	bool rowEmpty[capacity_];
	unsigned int rowMax[capacity_];
	unsigned int rowNextMax[capacity_];
	unsigned int rowNSuffixes[capacity_];
	unsigned int rowNSubnodes[capacity_];
	Node<DIM>* rowNode[capacity_];

	assert (!flushing_);
	flushing_ = true;
	const size_t currentIndex = nextIndex_;
	const size_t capacity = capacity_;
	const size_t n = min(currentIndex, capacity);
	assert (n > 0 && n <= capacity_);

	for (unsigned row = 0; row < n; ++row) {
		// spin until remaining insertions are done
		assert (insertCompleted_[row]);
		rowEmpty[row] = (getLcp(row, row) == -1u);
		rowNode[row] = NULL;
		rowMax[row] = -1u; // TODO not needed ?!
		rowNextMax[row] = -1u;
	}

	const unsigned int basicMsbIndex = WIDTH - suffixBits_ / DIM;
	bool hasMoreRows = true;
	while (hasMoreRows) {
		hasMoreRows = false;
		unsigned int maxRowMax = 0;

		// calculate maximum and number of occurrences per row
		// TODO no need to revisit rows that were not changed!
		for (unsigned row = 0; row < n; ++row) {
			if (rowEmpty[row]) continue;
			hasMoreRows = row != 0;

//			if (lastMaxRowMax == rowMax[row]) {
				// the row needs to be updated as the maximum changed

				rowMax[row] = -1u; // TODO possible similar to: (rowNextMax[row] == (-1u))? -1u : rowNextMax[row] - 1;
				rowNextMax[row] = -1u;

				for (unsigned column = 0; column < n; ++column) { // TODO how to save iterations?
					if (rowEmpty[column] || row == column) continue;

					const unsigned int currentLcp = getLcp(row, column);
					if ((1 + currentLcp) <= (1 + rowNextMax[row])) {
						// most common case: no need to update any values
						continue;
					} else if (currentLcp == rowMax[row] && rowNode[column]) {
						// found another LCP of the same length that has a subnode
						++rowNSubnodes[row];
					} else if (currentLcp == rowMax[row]) {
						// found another LCP of the same length without a subnode
						++rowNSuffixes[row];
					} else if ((1 + currentLcp) > (1 + rowMax[row])) {
						// found a higher LCP: store the new one
						rowNextMax[row] = rowMax[row];
						rowMax[row] = currentLcp;
						rowNSuffixes[row] = 0;
						rowNSubnodes[row] = 0;
						if (rowNode[row])	{++rowNSubnodes[row];}
						else 				{++rowNSuffixes[row];}
						if (rowNode[column]){++rowNSubnodes[row];}
						else 				{++rowNSuffixes[row];}
					} else { // TODO probably more likely and less expensive than case above
						assert ((1 + rowNextMax[row]) < (1 + currentLcp));
						// the second highest value needs to be updated
						rowNextMax[row] = currentLcp;
					}

					assert (rowNextMax[row] < rowMax[row] || rowNextMax[row] == (-1u));
				}
//			}

			assert ((1 + rowNextMax[row]) <= (1 + rowMax[row]));
//			assert (rowNSuffixes[row] + rowNSubnodes[row] <= (1uL << DIM));
			if ((1 + rowMax[row]) > (1 + maxRowMax)) {maxRowMax = rowMax[row];}
		}

#ifdef PRINT
		cout << "type\tmax\tnext\t#suff\t#node\t| i";
		for (unsigned row = 0; row < n; ++row) { cout << "\t"  << row;}
		cout << endl;
		// print the current LCP matrix
		for (unsigned row = 0; row < n; ++row) {
			if (rowNode[row]) { cout << "(node)\t"; } else { cout << "(suff)\t"; }
			if (rowMax[row] == -1u) { cout << "-1\t"; } else { cout << rowMax[row] << "\t"; }
			if (rowNextMax[row] == -1u) { cout << "-1\t"; } else {cout << rowNextMax[row] << "\t"; }
			cout << rowNSuffixes[row] << "\t" << rowNSubnodes[row] << "\t| " << row;
			for (unsigned col = 0; col < n; ++col) {
				cout << "\t";
				if (rowEmpty[row] || rowEmpty[col]) {cout << "-";}
				else if (getLcp(row, col) == (-1u)) {cout << "(-1)";}
				else {
					cout << getLcp(row, col);
					if (getLcp(row, col) == maxRowMax) {cout << "*";}
				}
			}
			cout << endl;
		}
		cout << endl;
#endif

		// create a new node for each row that was not ruled out yet
		const unsigned int index = basicMsbIndex + maxRowMax;
		assert (index < WIDTH);
		// current suffix bits: buffer suffix bits - current longest prefix bits - HC address bits (one node)
		const unsigned int suffixBits = suffixBits_ - DIM * (maxRowMax + 1);
		for (unsigned row = 0; row < (n - 1) && hasMoreRows; ++row) {
			if (rowEmpty[row] || maxRowMax != rowMax[row]) continue;

			setLcp(row, row, maxRowMax);
			assert ((1 + rowMax[row]) > (1 + rowNextMax[row]));
			// <nextMax>
			// <       current max        >
			// [ ------ |-  prefix -| DIM | .... ]
			const unsigned int prefixBits = (rowMax[row] - (1 + rowNextMax[row])) * DIM;
			const unsigned int nEntries = rowNSuffixes[row] + rowNSubnodes[row];
			assert (nEntries <= (1uL << DIM));
			Node<DIM>* currentNode = NodeTypeUtil<DIM>::
					template buildNodeWithSuffixes<WIDTH>(prefixBits, nEntries, rowNSuffixes[row], suffixBits);
			if (prefixBits > 0) {
				const unsigned int currentBits = suffixBits + DIM + prefixBits;
				assert (currentBits <= suffixBits_);
				MultiDimBitset<DIM>::duplicateBits(buffer_[row].values_,
						currentBits, prefixBits / DIM, currentNode->getPrefixStartBlock());
			}

			// insert all entries within this row into the new sub node
			for (unsigned column = row; column < n; ++column) {
				const unsigned int currentLcp = getLcp(row, column);
				if (rowEmpty[column] || currentLcp != rowMax[row]) continue;

				const unsigned long hcAddress = MultiDimBitset<DIM>::interleaveBits(buffer_[column].values_,index , DIM * WIDTH);
				assert (!(currentNode->lookup(hcAddress, true).exists));
				assert (currentNode->getNumberOfContents() < currentNode->getMaximumNumberOfContents());
				if (rowNode[column]) {
					// insert the subnode
					currentNode->insertAtAddress(hcAddress, rowNode[column]);
					 // TODO double???					setLcp(row, column, rowNextMax[row]);
				} else {
					// insert the suffix
					// TODO maybe sort before insertion so there is no need for swapping and sorting and shifting in LHC

					if (currentNode->canStoreSuffixInternally(suffixBits)) {
						unsigned long suffix = 0uL;
						MultiDimBitset<DIM>::removeHighestBits(buffer_[column].values_, DIM * WIDTH, index + 1, &suffix);
						currentNode->insertAtAddress(hcAddress, suffix, buffer_[column].id_);
					} else {
						assert (currentNode->canStoreSuffix(suffixBits) == 0);
						const pair<unsigned long*, unsigned int> suffixStartBlock = currentNode->reserveSuffixSpace(suffixBits);
						currentNode->insertAtAddress(hcAddress, suffixStartBlock.second, buffer_[column].id_);
						MultiDimBitset<DIM>::duplicateLowestBitsAligned(buffer_[column].values_, suffixBits, suffixStartBlock.first);
						assert(currentNode->lookup(hcAddress, true).suffixStartBlock == suffixStartBlock.first);
					}
				}

				rowEmpty[column] |= (row != column);
				rowNode[column] = currentNode;
			}

			setLcp(row, row, rowNextMax[row]); // TODO not needed?!
		}

		#ifndef NDEBUG
			// Validate round
			for (unsigned i = 0; i < n; ++i) {
				assert (!rowNode[i] || rowNode[i]->getNumberOfContents() >= 2);
			}
		#endif
	}

	assert ((this->node_->lookup(this->nodeHcAddress, true).exists)
			&& (this->node_->lookup(this->nodeHcAddress, true).hasSpecialPointer));
	this->node_->insertAtAddress(this->nodeHcAddress, rowNode[0]);

#ifndef NDEBUG
	// Validate all copied entries (except for the first one which was only copied partially)
	for (unsigned i = 0; i < n; ++i) {
		assert (insertCompleted_[i]);
		assert (i == 0 || rowEmpty[i]);
	}
#endif

	return rowNode[0];
}


#endif /* SRC_UTIL_ENTRYBUFFER_H_ */
