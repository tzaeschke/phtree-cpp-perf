/*
 * TNode.h
 *
 *  Created on: May 4, 2016
 *      Author: max
 */

#ifndef SRC_NODES_TNODE_H_
#define SRC_NODES_TNODE_H_

#include "../Entry.h"
#include "../iterators/NodeIterator.h"
#include "../nodes/NodeAddressContent.h"
#include "../util/MultiDimBitset.h"
#include "../nodes/Node.h"

template <unsigned int DIM>
class Visitor;
template <unsigned int DIM>
class SizeVisitor;
template <unsigned int DIM>
class PrefixSharingVisitor;
template <unsigned int DIM, unsigned int WIDTH>
class DynamicNodeOperationsUtil;

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class TNode : public Node<DIM> {
	// TODO remove friends and use getters and setters
	friend class NodeIterator<DIM>;
	friend class SizeVisitor<DIM>;
	friend class PrefixSharingVisitor<DIM>;
public:

	explicit TNode(size_t prefixLength);
	explicit TNode(TNode<DIM, PREF_BLOCKS>* other);
	virtual ~TNode() {}
	virtual std::ostream& output(std::ostream& os, size_t depth, size_t index, size_t totalBitLength) override;
	virtual NodeIterator<DIM>* begin() const = 0;
	virtual NodeIterator<DIM>* it(unsigned long hcAddress) const =0;
	virtual NodeIterator<DIM>* end() const = 0;
	virtual void accept(Visitor<DIM>* visitor, size_t depth, unsigned int index) override;
	virtual void recursiveDelete() = 0;
	// gets the number of contents: #suffixes + #subnodes
	virtual size_t getNumberOfContents() const = 0;
	virtual size_t getMaximumNumberOfContents() const = 0;
	virtual void lookup(unsigned long address, NodeAddressContent<DIM>& outContent, bool resolveSuffixIndex) const = 0;
	virtual void insertAtAddress(unsigned long hcAddress, uintptr_t pointer) =0;
	virtual void insertAtAddress(unsigned long hcAddress, unsigned int suffixStartBlockIndex, int id) = 0;
	virtual void insertAtAddress(unsigned long hcAddress, unsigned long startSuffixBlock, int id) = 0;
	virtual void insertAtAddress(unsigned long hcAddress, const Node<DIM>* const subnode) = 0;
	virtual Node<DIM>* adjustSize() = 0;

	size_t getMaxPrefixLength() const override;
	size_t getPrefixLength() const override;
	unsigned long* getPrefixStartBlock() override;
	const unsigned long* getFixPrefixStartBlock() const override;
	bool canStoreSuffixInternally(size_t nSuffixBits) const override;
	unsigned int canStoreSuffix(size_t nSuffixBits) const override;
	void setSuffixStorage(TSuffixStorage* suffixStorage) override;
	const TSuffixStorage* getSuffixStorage() const override;
	std::pair<unsigned long*, unsigned int> reserveSuffixSpace(size_t nSuffixBits) override;
	void freeSuffixSpace(size_t nSuffixBits, unsigned long* suffixStartBlock) override;
	void copySuffixStorageFrom(const Node<DIM>& other) override;

	unsigned long* getSuffixStartBlockPointerFromIndex(unsigned int index) const;
protected:
	size_t prefixBits_;
	TSuffixStorage* suffixes_;
	// TODO template block type
	unsigned long prefix_[PREF_BLOCKS];

	TSuffixStorage* getChangeableSuffixStorage() const override;
	virtual string getName() const =0;
};

#include <assert.h>
#include <stdexcept>
#include "../nodes/LHC.h"
#include "../util/SpatialSelectionOperationsUtil.h"
#include "../iterators/RangeQueryIterator.h"
#include "../iterators/NodeIterator.h"
#include "../nodes/TSuffixStorage.h"

using namespace std;

template <unsigned int DIM, unsigned int PREF_BLOCKS>
TNode<DIM, PREF_BLOCKS>::TNode(size_t prefixLength) : prefixBits_(prefixLength * DIM), suffixes_(NULL), prefix_() {
	assert (prefixBits_ <= PREF_BLOCKS * sizeof (unsigned long) * 8);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
TNode<DIM, PREF_BLOCKS>::TNode(TNode<DIM, PREF_BLOCKS>* other) : prefixBits_(other->prefixBits_), suffixes_(other->suffixes_) {
	assert (prefixBits_ <= PREF_BLOCKS * sizeof (unsigned long) * 8);
	for (unsigned i = 0; i < PREF_BLOCKS; ++i) {
		prefix_[i] = other->prefix_[i];
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
size_t TNode<DIM, PREF_BLOCKS>::getMaxPrefixLength() const {
	return PREF_BLOCKS * 8 * sizeof (unsigned long);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
size_t TNode<DIM, PREF_BLOCKS>::getPrefixLength() const {
	return prefixBits_ / DIM;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
unsigned long* TNode<DIM, PREF_BLOCKS>::getPrefixStartBlock() {
	return prefix_;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
const unsigned long* TNode<DIM, PREF_BLOCKS>::getFixPrefixStartBlock() const {
	return prefix_;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void TNode<DIM, PREF_BLOCKS>::accept(Visitor<DIM>* visitor, size_t depth, unsigned int index) {

	NodeIterator<DIM>* it;
	NodeIterator<DIM>* endIt = this->end();
	for (it = this->begin(); (*it) != *endIt; ++(*it)) {
		NodeAddressContent<DIM> content = *(*it);
		assert (content.exists);
		if (content.hasSubnode) {
			const size_t prefixLength = content.subnode->getPrefixLength();
			content.subnode->accept(visitor, depth + 1, index + prefixLength + 1);
		}
	}

	delete it;
	delete endIt;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
bool TNode<DIM, PREF_BLOCKS>::canStoreSuffixInternally(size_t nSuffixBits) const {
	// needs to store 2 bits for meta data and 32 bits for the ID
	return nSuffixBits <= 30;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
unsigned int TNode<DIM, PREF_BLOCKS>::canStoreSuffix(size_t nSuffixBits) const {
	assert(nSuffixBits > 0);

	if (suffixes_) {
		if (suffixes_->canStoreBits(nSuffixBits)) {
			return 0;
		} else {
			return suffixes_->getTotalBlocksToStoreAdditionalSuffix(nSuffixBits);
		}
	} else {
		const unsigned int suffixBlocks = 1 + (nSuffixBits - 1) / (8 * sizeof (unsigned long));
		return suffixBlocks;
	}
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void TNode<DIM, PREF_BLOCKS>::setSuffixStorage(TSuffixStorage* suffixes) {
	suffixes_ = suffixes;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
const TSuffixStorage* TNode<DIM, PREF_BLOCKS>::getSuffixStorage() const {
	return suffixes_;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
TSuffixStorage* TNode<DIM, PREF_BLOCKS>::getChangeableSuffixStorage() const {
	return suffixes_;
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
std::pair<unsigned long*, unsigned int> TNode<DIM, PREF_BLOCKS>::reserveSuffixSpace(size_t nSuffixBits) {
	assert (suffixes_ && suffixes_->canStoreBits(nSuffixBits));
	// assumes insertion order: reserve space -> fill suffix -> insert index
	assert (this->getNStoredSuffixes() == suffixes_->getNStoredSuffixes(nSuffixBits));
	return suffixes_->reserveBits(nSuffixBits);
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
unsigned long* TNode<DIM, PREF_BLOCKS>::getSuffixStartBlockPointerFromIndex(unsigned int index) const {
	assert (suffixes_);
	return suffixes_->getPointerFromIndex(index);
}

template<unsigned int DIM, unsigned int PREF_BLOCKS>
void TNode<DIM, PREF_BLOCKS>::freeSuffixSpace(size_t nSuffixBits,
		unsigned long* suffixStartBlock) {
	assert(suffixes_ && !suffixes_->empty());

	const unsigned int suffixStartBlockIndex = suffixes_->getIndexFromPointer(suffixStartBlock);
	const unsigned int lastStartBlockIndex = suffixes_->overrideBlocksWithLast(nSuffixBits, suffixStartBlockIndex);
	if (lastStartBlockIndex != 0) {
		NodeIterator<DIM>* it = this->begin();
		it->disableResolvingSuffixIndex();
		NodeIterator<DIM>* endIt = this->end();
		bool foundLastRef = false;
		for (; (*it != *endIt) && !foundLastRef; ++(*it)) {
			NodeAddressContent<DIM> content = (*(*it));
			if (content.exists && !content.hasSubnode
					&& content.suffixStartBlockIndex == lastStartBlockIndex) {
				assert(!content.directlyStoredSuffix);
				insertAtAddress(content.address, suffixStartBlockIndex, content.id);
				foundLastRef = true;
			}
		}

		assert(foundLastRef);
		delete it;
		delete endIt;
	} else {
		suffixes_->clearLast(nSuffixBits);
	}

	assert (this->getNStoredSuffixes() == suffixes_->getNStoredSuffixes(nSuffixBits));
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
void TNode<DIM, PREF_BLOCKS>::copySuffixStorageFrom(const Node<DIM>& other) {
	suffixes_ = other.getChangeableSuffixStorage();
}

template <unsigned int DIM, unsigned int PREF_BLOCKS>
ostream& TNode<DIM, PREF_BLOCKS>::output(std::ostream& os, size_t depth, size_t index, size_t totalBitLength) {
	os << this->getName() << " | prefix: ";
	MultiDimBitset<DIM>::output(os, prefix_, prefixBits_);
	os << endl;
	const size_t currentIndex = index + this->getPrefixLength() + 1;

	NodeIterator<DIM>* it;
	NodeIterator<DIM>* endIt = this->end();
	for (it = this->begin(); *it != *endIt; ++(*it)) {
		NodeAddressContent<DIM> content = (*(*it));
		for (size_t i = 0; i < depth; i++) {os << "-";}
		os << " " << content.address << ": ";

		if (content.hasSubnode) {
			// print subnode
			content.subnode->output(os, depth + 1, currentIndex, totalBitLength);
		} else {
			// print suffix
			os << " suffix: ";
			MultiDimBitset<DIM>::output(os, content.getSuffixStartBlock(), DIM * (totalBitLength - currentIndex));
			os << " (id: " << content.id << ")" << endl;
		}
	}

	delete it;
	delete endIt;
	return os;
}

#endif /* SRC_NODES_TNODE_H_ */
