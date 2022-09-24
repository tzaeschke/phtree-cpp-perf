/*
 * LHC.h
 *
 *  Created on: Feb 25, 2016
 *      Author: max
 */

#ifndef LHC_H_
#define LHC_H_

#include <map>
#include <vector>
#include <cstdint>
#include "../nodes/TNode.h"

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
class LHCIterator;

template <unsigned int DIM>
class AssertionVisitor;

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
class LHC: public TNode<DIM, PREF_BLOCKS> {
	friend class LHCIterator<DIM, PREF_BLOCKS, N>;
	friend class AssertionVisitor<DIM>;
	friend class SizeVisitor<DIM>;
public:
	explicit LHC(size_t prefixLength);
	virtual ~LHC();
	NodeIterator<DIM>* begin() const override;
	NodeIterator<DIM>* it(unsigned long hcAddress) const override;
	NodeIterator<DIM>* end() const override;
	void accept(Visitor<DIM>* visitor, size_t depth, unsigned int index) override;
	void recursiveDelete() override;
	size_t getNumberOfContents() const override;
	size_t getMaximumNumberOfContents() const override;
	void lookup(unsigned long address, NodeAddressContent<DIM>& outContent, bool resolveSuffixIndex) const override;
	void insertAtAddress(unsigned long hcAddress, uintptr_t pointer) override;
	void insertAtAddress(unsigned long hcAddress, unsigned int suffixStartBlockIndex, int id) override;
	void insertAtAddress(unsigned long hcAddress, unsigned long suffix, int id) override;
	void insertAtAddress(unsigned long hcAddress, const Node<DIM>* const subnode) override;
	Node<DIM>* adjustSize() override;

protected:
	string getName() const override;

private:
	static const unsigned long fullBlock = -1;
	static const unsigned int bitsPerBlock = sizeof (unsigned long) * 8;

	// map of valid HC addresses in ascending order
	// block : <-------------------- 64 --------------------><--- ...
	// bits  : <-   DIM    ->|<-   DIM    -> * N
	// N rows: [ hc address ] [ hc address ] ...
	unsigned long addresses_[1 + ((N * DIM) - 1) / bitsPerBlock];
	// stores flags in 2 lowest bits per reference:
	// meaning of flags: isPointer | isSuffix
	// 00 - special pointer
	// 01 - the entry directly stores a suffix and the ID
	// 10 - the entry holds a reference to a subnode
	// 11 - the entry holds the index of the suffix and the ID
	std::uintptr_t references_[N];
	// number of actually filled rows: 0 <= m <= N
	unsigned int m;

	// <found?, index, hasSub?>
	void lookupAddress(unsigned long hcAddress, bool* outExists, unsigned int* outIndex) const;
	void lookupIndex(unsigned int index, unsigned long* outHcAddress) const;
	void fillLookupContent(NodeAddressContent<DIM>& outContent, uintptr_t reference, bool resolveSuffixIndex) const;
	inline void addRow(unsigned int index, unsigned long hcAddress, std::uintptr_t reference);
	inline void insertAddress(unsigned int index, unsigned long hcAddress);
	inline void interpretReference(std::uintptr_t ref, bool* isPointer, bool* isSuffix) const;
};

#include <assert.h>
#include <utility>
#include "../nodes/AHC.h"
#include "../nodes/LHC.h"
#include "../iterators/LHCIterator.h"
#include "../visitors/Visitor.h"
#include "../util/NodeTypeUtil.h"

using namespace std;
template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
LHC<DIM, PREF_BLOCKS, N>::LHC(size_t prefixLength) : TNode<DIM, PREF_BLOCKS>(prefixLength),
	addresses_(), references_(), m(0) {
	assert (N > 0 && m >= 0 && m <= N);
	assert (N <= (1 << DIM));
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
LHC<DIM, PREF_BLOCKS, N>::~LHC() {
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::recursiveDelete() {
	bool isPointer, isSuffix;
	const unsigned long flagMask = ~(3uL);
	for (unsigned int i = 0; i < m; ++i) {
		interpretReference(references_[i], &isPointer, &isSuffix);
		if (isPointer && !isSuffix) {
			Node<DIM>* subnode = reinterpret_cast<Node<DIM>*>(references_[i] & flagMask);
			subnode->recursiveDelete();
		}
	}

	if (this->suffixes_) { delete this->suffixes_; }
	delete this;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
string LHC<DIM,PREF_BLOCKS, N>::getName() const {
	return "LHC";
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM,PREF_BLOCKS, N>::lookupIndex(unsigned int index, unsigned long* outHcAddress) const {
	assert (index < m && m <= N);
	assert (DIM <= bitsPerBlock);

	const unsigned int firstBit = index * DIM;
	const unsigned int lastBit = (index + 1) * DIM;
	const unsigned int firstBitBlockIndex = firstBit / bitsPerBlock;
	const unsigned int lastBitBlockIndex = lastBit / bitsPerBlock;
	const unsigned int firstBitIndex = firstBit % bitsPerBlock;
	const unsigned int lastBitIndex = lastBit % bitsPerBlock;

	const unsigned long firstBlock = addresses_[firstBitBlockIndex];
	assert (lastBitBlockIndex - firstBitBlockIndex <= 1);
	if (firstBitBlockIndex == lastBitBlockIndex || lastBitIndex == 0) {
		// all required bits are in one block
		const unsigned long singleBlockAddressMask = (1uL << DIM) - 1uL;
		const unsigned long extracted = firstBlock >> firstBitIndex;
		(*outHcAddress) = extracted & singleBlockAddressMask;
	} else {
		// the address is split into two blocks
		const unsigned long secondBlock = addresses_[lastBitBlockIndex];
		assert (0 < lastBitIndex && lastBitIndex < DIM);
		const unsigned int firstBlockBits = bitsPerBlock - firstBitIndex;
		const unsigned int secondBlockBits = DIM - firstBlockBits;
		const unsigned long secondBlockMask = (1uL << secondBlockBits) - 1uL;
		(*outHcAddress) = (firstBlock >> firstBitIndex)
						| ((secondBlock & secondBlockMask) << firstBlockBits);
	}

	assert (*outHcAddress < (1uL << DIM));
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM,PREF_BLOCKS, N>::insertAddress(unsigned int index, unsigned long hcAddress) {
	assert (index < m && m <= N);
	assert (DIM <= bitsPerBlock);

	const unsigned int firstBit = index * DIM;
	const unsigned int lastBit = (index + 1) * DIM;
	const unsigned int firstBlockIndex = firstBit / bitsPerBlock;
	const unsigned int secondBlockIndex = lastBit / bitsPerBlock;
	const unsigned int firstBitIndex = firstBit % bitsPerBlock;
	const unsigned int lastBitIndex = lastBit % bitsPerBlock;
	assert (secondBlockIndex - firstBlockIndex <= 1);
	assert (lastBit - firstBit == DIM);

	const unsigned long firstBlockAddressMask = (1uL << DIM) - 1;
	// [   (   x   )   ]
	addresses_[firstBlockIndex] &= ~(firstBlockAddressMask << firstBitIndex);
	addresses_[firstBlockIndex] |= hcAddress << firstBitIndex;

	if (firstBlockIndex != secondBlockIndex && lastBitIndex != 0) {
		assert (secondBlockIndex - firstBlockIndex == 1);
		// the required bits are in two consecutive blocks
		//       <---------->
		// [     ( x  ] [ y )    ]
		assert (DIM > lastBitIndex);
		const unsigned int firstBlockBits = bitsPerBlock - firstBitIndex;
		const unsigned int secondBlockBits = DIM - firstBlockBits;
		assert (0 < secondBlockBits && secondBlockBits < DIM);
		const unsigned long secondBlockMask = fullBlock << secondBlockBits;
		const unsigned long remainingHcAddress = hcAddress >> firstBlockBits;
		assert (remainingHcAddress < (1uL << secondBlockBits));
		addresses_[secondBlockIndex] &= secondBlockMask;
		addresses_[secondBlockIndex] |= remainingHcAddress;
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::lookupAddress(unsigned long hcAddress, bool* outExists,
		unsigned int* outIndex) const {
	if (m == 0) {
		(*outExists) = false;
		(*outIndex) = 0;
		return;
	}

	// perform binary search in range [0, m)
	unsigned int l = 0;
	unsigned int r = m;
	unsigned long currentHcAddress = -1;
	while (l < r) {
		// check interval [l, r)
		const unsigned int middle = (l + r) / 2;
		assert (0 <= middle && middle < m);
		lookupIndex(middle, &currentHcAddress);
		if (currentHcAddress < hcAddress) {
			l = middle + 1;
		} else if (currentHcAddress > hcAddress) {
			r = middle;
		} else {
			// found the correct index
			(*outExists) = true;
			(*outIndex) = middle;
			return;
		}
	}

	// did not find the entry so set the position it should have
	assert (l - r == 0);
	(*outExists) = false;
	(*outIndex) = r;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::interpretReference(std::uintptr_t ref, bool* isPointer, bool* isSuffix) const {
	(*isSuffix) = ref & 1uL;
	(*isPointer) = (ref >> 1uL) & 1uL;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::fillLookupContent(NodeAddressContent<DIM>& outContent,
		uintptr_t reference, bool resolveSuffixIndex) const {
	assert (outContent.exists);
	bool isPointer, isSuffix;
	interpretReference(reference, &isPointer, &isSuffix);
	outContent.hasSubnode = isPointer && !isSuffix;
	outContent.directlyStoredSuffix = !isPointer && isSuffix;
	outContent.hasSpecialPointer = !isPointer && !isSuffix;

	if (outContent.exists) {
		const unsigned long flagMask = ~(3uL);
		if (outContent.hasSubnode) {
			outContent.subnode = reinterpret_cast<Node<DIM>*>(reference & flagMask);
		} else if (outContent.hasSpecialPointer) {
			assert ((reference & flagMask) == reference);
			outContent.specialPointer = reference;
		} else {
			const unsigned long suffixAndId = reinterpret_cast<unsigned long>(reference);
			const unsigned long suffixMask = (-1uL) >> 32;
			const unsigned int suffixPart = (suffixAndId & suffixMask) >> 2;
			if (outContent.directlyStoredSuffix) {
				outContent.suffix = suffixPart;
			} else {
				if (resolveSuffixIndex) {
					outContent.suffixStartBlock = this->getSuffixStartBlockPointerFromIndex(suffixPart);
				} else {
					outContent.suffixStartBlockIndex = suffixPart;
				}
			}

			const unsigned long idPart = (suffixAndId & (~suffixMask)) >> 32;
			assert (idPart < (1uL << 32));
			outContent.id = idPart;
		}
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::lookup(unsigned long address, NodeAddressContent<DIM>& outContent, bool resolveSuffixIndex) const {
	assert (address < 1uL << DIM);
	outContent.address = address;
	unsigned int index = m;
	lookupAddress(address, &outContent.exists, &index);
	if (outContent.exists) {
		fillLookupContent(outContent, references_[index], resolveSuffixIndex);
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::addRow(unsigned int index, unsigned long newHcAddress,
		uintptr_t newReference) {

	assert (!((NodeAddressContent<DIM>)Node<DIM>::lookup(newHcAddress, true)).exists);

	++m;
	if (index != m - 1) {
		assert (index < m && m <= N);
		// move all contents as of the given index
		// (copying from lower index because of forward prefetching)
		uintptr_t lastRef = references_[index];
		unsigned long lastAddress = 0;
		lookupIndex(index, &lastAddress);
		assert (lastAddress > newHcAddress);
		unsigned long tmpAddress = 0;

		for (unsigned i = index + 1; i < m - 1; ++i) {
			const uintptr_t tmpRef = references_[i];
			lookupIndex(i, &tmpAddress);
			assert (tmpAddress > newHcAddress);
			references_[i] = lastRef;
			insertAddress(i, lastAddress);
			lastRef = tmpRef;
			lastAddress = tmpAddress;
		}

		references_[m - 1] = lastRef;
		insertAddress(m - 1, lastAddress);
	}

	// insert the new entry at the freed index
	references_[index] = newReference;
	insertAddress(index, newHcAddress);
	assert (m <= N);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::insertAtAddress(unsigned long hcAddress, uintptr_t pointer) {
	assert (hcAddress < 1uL << DIM);
	assert ((pointer & 3) == 0);
	// format: [ pointer (62) | flags - 00 (2) ]

	unsigned int index = m;
	bool exists;
	lookupAddress(hcAddress, &exists, &index);
	assert (index <= m);

	if (exists) {
		// replace the contents at the address
		references_[index] = pointer;
	} else {
		// add a new entry
		assert (m < N && "the maximum number of entries must not have been reached");
		addRow(index, hcAddress, pointer);
	}

	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).hasSpecialPointer);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).specialPointer == pointer);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::insertAtAddress(unsigned long hcAddress, unsigned int suffixStartBlockIndex, int id) {
	assert (hcAddress < 1uL << DIM);
	assert (suffixStartBlockIndex < (1uL << 30));
	// format: [ ID (32) | suffix index (30) | flags - 11 (2) ]

	unsigned int index = m;
	bool exists;
	lookupAddress(hcAddress, &exists, &index);
	assert (index <= m);
	const unsigned long upperId = id;
	const unsigned long suffixStartBlockIndexExtended = (upperId << 32) | (suffixStartBlockIndex << 2) | 3;
	const uintptr_t reference = reinterpret_cast<uintptr_t>(suffixStartBlockIndexExtended);

	if (exists) {
		// replace the contents at the address
		references_[index] = reference;
	} else {
		// add a new entry
		assert (m < N && "the maximum number of entries must not have been reached");
		addRow(index, hcAddress, reference);
	}

	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).id == id);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (!((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).hasSubnode);
	assert (!((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).directlyStoredSuffix);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::insertAtAddress(unsigned long hcAddress, unsigned long suffix, int id) {
	assert (hcAddress < 1uL << DIM);
	assert (suffix < (1uL << 30));
	// format: [ ID (32) | suffix (30) | flags - 01 (2) ]

	unsigned int index = m;
	bool exists;
	lookupAddress(hcAddress, &exists, &index);
	assert (index <= m);
	const unsigned long upperId = id;
	const uintptr_t reference = reinterpret_cast<uintptr_t>((upperId << 32) | (suffix << 2) | 1);

	if (exists) {
		// replace the contents at the address
		references_[index] = reference;
	} else {
		// add a new entry
		assert (m < N && "the maximum number of entries must not have been reached");
		addRow(index, hcAddress, reference);
	}

	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).id == id);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).directlyStoredSuffix);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).suffix == suffix);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::insertAtAddress(unsigned long hcAddress, const Node<DIM>* const subnode) {
	assert (subnode);
	assert (hcAddress < 1uL << DIM);
	// format: [ subnode reference (62) | flags - 10 (2) ]

	unsigned int index = m;
	bool exists;
	lookupAddress(hcAddress, &exists, &index);
	assert (index <= m);
	const uintptr_t subRef = reinterpret_cast<uintptr_t>(subnode);
	assert ((subRef & 3) == 0);
	const uintptr_t reference = reinterpret_cast<uintptr_t>(subRef | 2);

	if (exists) {
		// replace the contents at the address
		references_[index] = reference;
	} else {
		// add a new entry
		assert (m < N && "the maximum number of entries must not have been reached");
		addRow(index, hcAddress, reference);
	}

#ifndef NDEBUG
	NodeAddressContent<DIM> content = Node<DIM>::lookup(hcAddress, true);
	assert (content.address == hcAddress);
	assert (content.hasSubnode);
	assert (content.subnode == subnode);
#endif
}


template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
Node<DIM>* LHC<DIM, PREF_BLOCKS, N>::adjustSize() {
	// TODO put this method into insert method instead
	if (m <= N) {
		return this;
	} else {
		return NodeTypeUtil<DIM>::copyIntoLargerNode(N + 1, this);
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
NodeIterator<DIM>* LHC<DIM, PREF_BLOCKS, N>::begin() const {
	LHCIterator<DIM, PREF_BLOCKS, N>* it = new LHCIterator<DIM, PREF_BLOCKS, N>(*this);
	if (m == 0) it->setToEnd();
	else it->setToBegin();
	return it;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
NodeIterator<DIM>* LHC<DIM, PREF_BLOCKS, N>::it(unsigned long hcAddress) const {
	return new LHCIterator<DIM, PREF_BLOCKS, N>(hcAddress, *this);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
NodeIterator<DIM>* LHC<DIM, PREF_BLOCKS, N>::end() const {
	NodeIterator<DIM>* it = new LHCIterator<DIM, PREF_BLOCKS, N>(*this);
	it->setToEnd();
	return it;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
size_t LHC<DIM, PREF_BLOCKS, N>::getNumberOfContents() const {
	return m;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
size_t LHC<DIM, PREF_BLOCKS, N>::getMaximumNumberOfContents() const {
	return N;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
void LHC<DIM, PREF_BLOCKS, N>::accept(Visitor<DIM>* visitor, size_t depth, unsigned int index) {
	visitor->visit(this, depth, index);
	TNode<DIM, PREF_BLOCKS>::accept(visitor, depth, index);
}

#endif /* LHC_H_ */
