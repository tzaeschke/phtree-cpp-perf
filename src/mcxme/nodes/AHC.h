/*
 * AHC.h
 *
 *  Created on: Feb 25, 2016
 *      Author: max
 */

#ifndef SRC_AHC_H_
#define SRC_AHC_H_

#include <vector>
#include <iostream>
#include <cstdint>
#include "../nodes/TNode.h"
#include "../iterators/NodeIterator.h"
#include "../nodes/NodeAddressContent.h"

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class AHCIterator;

template <unsigned int DIM>
class AssertionVisitor;

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class AHC: public TNode<DIM, PREF_BLOCKS> {
	friend class AHCIterator<DIM, PREF_BLOCKS>;
	friend class AssertionVisitor<DIM>;
	friend class SizeVisitor<DIM>;
public:
	explicit AHC(size_t prefixLength);
	explicit AHC(TNode<DIM, PREF_BLOCKS>* node);
	virtual ~AHC();
	NodeIterator<DIM>* begin() const override;
	NodeIterator<DIM>* it(unsigned long hcAddress) const override;
	NodeIterator<DIM>* end() const override;
	void accept(Visitor<DIM>* visitor, size_t depth, unsigned int index) override;
	void recursiveDelete() override;
	size_t getNumberOfContents() const override;
	size_t getMaximumNumberOfContents() const override;
	void lookup(unsigned long address, NodeAddressContent<DIM>& outContent, bool resolveSuffixIndex) const override;
	void insertAtAddress(unsigned long hcAddress, uintptr_t pointer) override;
	void insertAtAddress(unsigned long hcAddress, unsigned long suffix, int id) override;
	void insertAtAddress(unsigned long hcAddress, unsigned int suffixStartBlockIndex, int id) override;
	void insertAtAddress(unsigned long hcAddress, const Node<DIM>* const subnode) override;
	Node<DIM>* adjustSize() override;

protected:
	string getName() const override;

private:

	size_t nContents;

	// stores flags in 2 lowest bits per reference: isPointer | isSuffix
	// 00 & reference == 0	- entry does not exist
	// 00 & reference != 0	- entry is a special pointer
	// 01 					- the entry directly stores a suffix and the ID
	// 10 					- the entry holds a reference to a subnode
	// 11 					- the entry holds the index of the suffix and the ID
	std::uintptr_t references_[1 << DIM];
	static const unsigned long refMask = (-1uL) << 2uL; // mask to remove the 2 flag bits

	void inline getRef(unsigned long hcAddress, bool* exists, bool* hasSub,
			bool* directlyStoredSuffix, bool* isSpecial, std::uintptr_t* ref) const;
};

#include <assert.h>
#include "../nodes/AHC.h"
#include "../iterators/AHCIterator.h"
#include "../iterators/NodeIterator.h"
#include "../nodes/NodeAddressContent.h"
#include "../visitors/Visitor.h"

using namespace std;

template <unsigned int DIM, unsigned int PREF_BLOCKS>
AHC<DIM, PREF_BLOCKS>::AHC(size_t prefixLength) : TNode<DIM, PREF_BLOCKS>(prefixLength), nContents(0), references_() {
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
AHC<DIM, PREF_BLOCKS>::AHC(TNode<DIM, PREF_BLOCKS>* other) : TNode<DIM, PREF_BLOCKS>(other), nContents(0), references_() {

	// TODO use more efficient way to convert LHC->AHC
	NodeIterator<DIM>* it;
	NodeIterator<DIM>* endIt = other->end();
	for (it = other->begin(); (*it) != *endIt; ++(*it)) {
		NodeAddressContent<DIM> content = *(*it);
		assert (content.exists);
		if (content.hasSubnode) {
			insertAtAddress(content.address, content.subnode);
		} else if (content.directlyStoredSuffix) {
			assert (content.suffix < sizeof (int) * 8 - 2);
			insertAtAddress(content.address, content.suffix, content.id);
		} else {
			assert (content.suffixStartBlock);
			insertAtAddress(content.address, content.suffixStartBlock, content.id);
		}
	}

	delete it;
	delete endIt;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
AHC<DIM, PREF_BLOCKS>::~AHC() {
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::recursiveDelete() {
	bool hasSubnode = false;
	bool filled = false;
	bool isDirectlyStoredSuffix = false;
	bool isSpecial = false;
	uintptr_t ref = 0;
	for (size_t i = 0; i < 1uL << DIM; ++i) {
		getRef(i, &filled, &hasSubnode, &isDirectlyStoredSuffix, &isSpecial, &ref);
		assert (!isSpecial);
		if (filled && hasSubnode) {
			Node<DIM>* subnode = reinterpret_cast<Node<DIM>*>(ref);
			assert (subnode);
			subnode->recursiveDelete();
		}
	}

	if (this->suffixes_) { delete this->suffixes_; }
	delete this;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
string AHC<DIM,PREF_BLOCKS>::getName() const {
	return "AHC";
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM,PREF_BLOCKS>::getRef(unsigned long hcAddress, bool* exists, bool* hasSub,
		bool* directlyStoredSuffix, bool* isSpecial, std::uintptr_t* ref) const {
	assert (hcAddress < 1 << DIM);
	assert (exists && hasSub && ref);

	const uintptr_t reference = references_[hcAddress];
	const bool isSuffix = reference & 1;
	const bool isPointer = (reference >> 1) & 1;

	(*ref) = reference & refMask;
	(*exists) = reference != 0;
	assert ((*exists) || reference == 0);
	(*hasSub) = isPointer && !isSuffix;
	(*directlyStoredSuffix) = !isPointer && isSuffix;
	(*isSpecial) = !isSuffix && !isPointer && (reference != 0);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
size_t AHC<DIM, PREF_BLOCKS>::getNumberOfContents() const {
	return nContents;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
size_t AHC<DIM, PREF_BLOCKS>::getMaximumNumberOfContents() const {
	return 1uL << DIM;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::lookup(unsigned long address, NodeAddressContent<DIM>& outContent, bool resolveSuffixIndex) const {
	assert (address < 1uL << DIM);

	uintptr_t ref;
	outContent.address = address;
	getRef(address, &outContent.exists, &outContent.hasSubnode,
			&outContent.directlyStoredSuffix, &outContent.hasSpecialPointer,
			&ref);

	if (outContent.exists) {
		if (outContent.hasSubnode) {
			outContent.subnode = reinterpret_cast<Node<DIM>*>(ref);
		} else if (outContent.hasSpecialPointer) {
			outContent.specialPointer = ref;
		} else {
			const unsigned long suffixAndId = reinterpret_cast<unsigned long>(ref);
			const unsigned long suffixMask = (-1uL) >> 32;
			// convert the stored suffix and shift back where meta flags were stored
			const unsigned int suffixPart = (suffixAndId & suffixMask) >> 2;
			const unsigned long idPart = (suffixAndId & (~suffixMask)) >> 32;
			assert (idPart < (1uL << 32));
			outContent.id = idPart;

			if (outContent.directlyStoredSuffix) {
				outContent.suffix = suffixPart;
			} else {
				if (resolveSuffixIndex) {
					outContent.suffixStartBlock = this->getSuffixStartBlockPointerFromIndex(suffixPart);
				} else {
					outContent.suffixStartBlockIndex = suffixPart;
				}
			}
		}
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::insertAtAddress(unsigned long hcAddress, unsigned long  suffix, int id) {
	assert (hcAddress < 1ul << DIM);
	assert (sizeof (unsigned long) == sizeof (uintptr_t));
	assert (suffix < (1uL << (sizeof(int) * 8uL - 2uL)));

	bool exists;
	bool hasSubnode;
	bool isDirectlyStoredSuffix;
	bool isSpecial;
	uintptr_t ref;
	getRef(hcAddress, &exists, &hasSubnode, &isDirectlyStoredSuffix, &isSpecial, &ref);

	if (!exists) {
		nContents++;
		assert ((references_[hcAddress] & 3) == 0);
	}

	// layout: [ ID (32) | suffix bits (30) | flags (2)]
	// 01 -> internally stored suffix
	// need to shift the suffix in order to have enough space for the meta flags
	const unsigned long extendedId = id;
	const unsigned long suffixShifted = 1 | (suffix << 2) | (extendedId << 32);
	references_[hcAddress] = reinterpret_cast<uintptr_t>(suffixShifted);


	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).id == id);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).directlyStoredSuffix);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).suffix == suffix);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::insertAtAddress(unsigned long hcAddress, unsigned int suffixStartBlockIndex, int id) {
	assert (hcAddress < 1ul << DIM);
	assert (suffixStartBlockIndex < (1uL << (8 * sizeof (int) - 2)));

	bool exists;
	bool hasSubnode;
	bool isDirectlyStoredSuffix;
	bool isSpecial;
	uintptr_t ref;
	getRef(hcAddress, &exists, &hasSubnode, &isDirectlyStoredSuffix, &isSpecial, &ref);

	if (!exists) {
		nContents++;
		assert ((references_[hcAddress] & 3) == 0);
	}

	const unsigned long suffixStartBlockIndexExtended = suffixStartBlockIndex;
	const unsigned long extendedId = id;
	const uintptr_t storedRef = reinterpret_cast<uintptr_t>(suffixStartBlockIndexExtended << 2);
	assert (((storedRef & 3) == 0) && "last 2 bits must not be set!");
	// layout: [ ID (32) | suffix index (30) | flags (2)]
	// 11 -> pointer to suffix
	references_[hcAddress] = 3 | storedRef | (extendedId << 32);


	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).id == id);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (!((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).hasSubnode);
	assert (!((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).directlyStoredSuffix);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::insertAtAddress(unsigned long hcAddress, const Node<DIM>* const subnode) {
	assert (hcAddress < 1ul << DIM);

	bool exists;
	bool hasSubnode;
	bool isDirectlyStoredSuffix;
	bool isSpecial;
	uintptr_t ref;
	getRef(hcAddress, &exists, &hasSubnode, &isDirectlyStoredSuffix, &isSpecial, &ref);

	if (!exists) {
		nContents++;
		assert ((references_[hcAddress] & 3) == 0);
	}

	assert (((reinterpret_cast<uintptr_t>(subnode) & 3) == 0) && "last 2 bits must not be set!");
	// 10 -> pointer to subnode
	references_[hcAddress] = 2 | reinterpret_cast<uintptr_t>(subnode);

	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).hasSubnode);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).subnode == subnode);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::insertAtAddress(unsigned long hcAddress, uintptr_t pointer) {
	assert (hcAddress < 1ul << DIM);
	assert ((pointer & 3) == 0 && "last 2 bits must not be set!");
	assert (pointer != 0 && "cannot store null pointers");

	bool exists;
	bool hasSubnode;
	bool isDirectlyStoredSuffix;
	bool isSpecial;
	uintptr_t ref;
	getRef(hcAddress, &exists, &hasSubnode, &isDirectlyStoredSuffix, &isSpecial, &ref);

	if (!exists) {
		nContents++;
	}

	// 00 & pointer != 0 -> special pointer
	references_[hcAddress] = pointer;

	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).address == hcAddress);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).hasSpecialPointer);
	assert (((NodeAddressContent<DIM>)Node<DIM>::lookup(hcAddress, true)).specialPointer == pointer);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
Node<DIM>* AHC<DIM, PREF_BLOCKS>::adjustSize() {
	// TODO currently there is no need to switch from AHC to LHC because there is no delete function
	return this;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
NodeIterator<DIM>* AHC<DIM, PREF_BLOCKS>::begin() const {
	AHCIterator<DIM, PREF_BLOCKS>* it = new AHCIterator<DIM, PREF_BLOCKS>(*this);
	it->setToBegin();
	return it;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
NodeIterator<DIM>* AHC<DIM, PREF_BLOCKS>::it(unsigned long hcAddress) const {
	return new AHCIterator<DIM, PREF_BLOCKS>(hcAddress, *this);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
NodeIterator<DIM>* AHC<DIM, PREF_BLOCKS>::end() const {
	NodeIterator<DIM>* it = new AHCIterator<DIM, PREF_BLOCKS>(*this);
	it->setToEnd();
	return it;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void AHC<DIM, PREF_BLOCKS>::accept(Visitor<DIM>* visitor, size_t depth, unsigned int index)  {
	visitor->visit(this, depth, index);
	TNode<DIM, PREF_BLOCKS>::accept(visitor, depth, index);
}

#endif /* SRC_AHC_H_ */
