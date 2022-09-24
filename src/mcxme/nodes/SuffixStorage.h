/*
 * SuffixStorage.h
 *
 *  Created on: Jun 26, 2016
 *      Author: max
 */

#ifndef SRC_NODES_SUFFIXSTORAGE_H_
#define SRC_NODES_SUFFIXSTORAGE_H_


template <unsigned int DIM>
class SizeVisitor;

template <unsigned int SUFFIX_BLOCKS>
class SuffixStorage : public TSuffixStorage {
	template <unsigned int D>
	friend class SizeVisitor;
public:
	SuffixStorage();
	virtual ~SuffixStorage() {};
	bool canStoreBits(size_t nBitsToStore) const override;
	unsigned int getTotalBlocksToStoreAdditionalSuffix(size_t nSuffixBits) const override;
	std::pair<unsigned long*, unsigned int> reserveBits(size_t nBits) override;
	unsigned int overrideBlocksWithLast(size_t nBits, unsigned int overrideStartBlockIndex) override;
	void copyFrom(const TSuffixStorage& other) override;
	void clear() override;
	void clearLast(size_t nBits) override;
	unsigned int getNMaxStorageBlocks() const override;
	unsigned int getNCurrentStorageBlocks() const override;
	unsigned long getBlock(unsigned int index) const override;
	bool empty() const override;
	unsigned long* getPointerFromIndex(unsigned int index) const override;
	unsigned int getIndexFromPointer(unsigned long* pointer) const override;
	size_t getByteSize() const override;

private:
	unsigned int currentBlock;
//	unsigned long occupiedBlocks[1 + (SUFFIX_BLOCKS - 1) / (8 * sizeof (unsigned long))];
	unsigned long suffixBlocks[SUFFIX_BLOCKS];
};

using namespace std;
#include <assert.h>

template <unsigned int SUFFIX_BLOCKS>
SuffixStorage<SUFFIX_BLOCKS>::SuffixStorage() : currentBlock(0), suffixBlocks() { }

template <unsigned int SUFFIX_BLOCKS>
unsigned int  SuffixStorage<SUFFIX_BLOCKS>::getNMaxStorageBlocks() const {
	return SUFFIX_BLOCKS;
}

template <unsigned int SUFFIX_BLOCKS>
unsigned int SuffixStorage<SUFFIX_BLOCKS>::getNCurrentStorageBlocks() const {
	return currentBlock;
}

template <unsigned int SUFFIX_BLOCKS>
bool SuffixStorage<SUFFIX_BLOCKS>::empty() const {
	return currentBlock == 0;
}

template <unsigned int SUFFIX_BLOCKS>
unsigned long SuffixStorage<SUFFIX_BLOCKS>::getBlock(unsigned int index) const {
	return suffixBlocks[index];
}

template <unsigned int SUFFIX_BLOCKS>
unsigned long* SuffixStorage<SUFFIX_BLOCKS>::getPointerFromIndex(unsigned int index) const {
	return const_cast<unsigned long*>(suffixBlocks + index);
}

template <unsigned int SUFFIX_BLOCKS>
unsigned int SuffixStorage<SUFFIX_BLOCKS>::getIndexFromPointer(unsigned long* pointer) const {
	const size_t pointerRepr = (size_t)pointer;
	const size_t startRepr = (size_t)suffixBlocks;
	const size_t diff = pointerRepr - startRepr;
	assert (diff < (1uL << (8 * sizeof (unsigned int))));
	const unsigned int index = diff / (sizeof (unsigned long));
	assert (index < currentBlock);
	assert ((*pointer) == suffixBlocks[index]);
	return index;
}

template <unsigned int SUFFIX_BLOCKS>
void SuffixStorage<SUFFIX_BLOCKS>::copyFrom(const TSuffixStorage& other) {
	assert (SUFFIX_BLOCKS >= other.getNCurrentStorageBlocks());
	currentBlock = other.getNCurrentStorageBlocks();
	for (unsigned int block = 0; block < currentBlock; ++block) {
		suffixBlocks[block] = other.getBlock(block);
	}
}

template <unsigned int SUFFIX_BLOCKS>
unsigned int SuffixStorage<SUFFIX_BLOCKS>::overrideBlocksWithLast(size_t nBits, unsigned int overrideStartBlockIndex) {
	assert (overrideStartBlockIndex < currentBlock);
	const size_t nBlocks = 1 + (nBits - 1) / (8 * sizeof (unsigned long));
	assert (overrideStartBlockIndex % nBlocks == 0);
	const unsigned int lastStartBlockIndex = currentBlock - nBlocks;
	if (lastStartBlockIndex == 0 || overrideStartBlockIndex == lastStartBlockIndex) {
		return 0;
	}

	assert (overrideStartBlockIndex < lastStartBlockIndex && lastStartBlockIndex < currentBlock);
	for (unsigned block = 0; block < nBlocks; ++block) {
		suffixBlocks[overrideStartBlockIndex + block] = suffixBlocks[lastStartBlockIndex + block];
		suffixBlocks[lastStartBlockIndex + block] = 0;
	}

	currentBlock -= nBlocks;
	return lastStartBlockIndex;
}

template <unsigned int SUFFIX_BLOCKS>
void SuffixStorage<SUFFIX_BLOCKS>::clear() {
	for (unsigned i = 0; i < currentBlock; ++i) {
		suffixBlocks[i] = 0;
	}

	currentBlock = 0;
}

template <unsigned int SUFFIX_BLOCKS>
void SuffixStorage<SUFFIX_BLOCKS>::clearLast(size_t nBits) {
	const size_t nBlocks = 1 + (nBits - 1) / (8 * sizeof (unsigned long));
	const unsigned int lastStartBlockIndex = currentBlock - nBlocks;
	for (unsigned i = lastStartBlockIndex; i < currentBlock; ++i) {
		suffixBlocks[i] = 0;
	}

	currentBlock -= nBlocks;
}

template <unsigned int SUFFIX_BLOCKS>
unsigned int SuffixStorage<SUFFIX_BLOCKS>::getTotalBlocksToStoreAdditionalSuffix(size_t nSuffixBits) const {
	const size_t nBlocks = 1 + (nSuffixBits - 1) / (8 * sizeof (unsigned long));
	assert (SUFFIX_BLOCKS < currentBlock + nBlocks);
	return currentBlock + nBlocks;
}

template <unsigned int SUFFIX_BLOCKS>
bool SuffixStorage<SUFFIX_BLOCKS>::canStoreBits(size_t nBitsToStore) const {
	const size_t nBlocks = 1 + (nBitsToStore - 1) / (8 * sizeof (unsigned long));
	return SUFFIX_BLOCKS >=  currentBlock + nBlocks;
}

template <unsigned int SUFFIX_BLOCKS>
pair<unsigned long*, unsigned int> SuffixStorage<SUFFIX_BLOCKS>::reserveBits(size_t nBits) {
	assert (canStoreBits(nBits));
	const size_t nBlocks = 1 + (nBits - 1) / (8 * sizeof (unsigned long));
	unsigned long* reservedStartBlock = suffixBlocks + currentBlock;
	unsigned int startBlockIndex = currentBlock;
	currentBlock += nBlocks;
	return pair<unsigned long*, unsigned int>(reservedStartBlock, startBlockIndex);
}

template <unsigned int SUFFIX_BLOCKS>
size_t SuffixStorage<SUFFIX_BLOCKS>::getByteSize() const {
	size_t byteSize = sizeof (currentBlock);
	byteSize += sizeof (suffixBlocks);
	return byteSize;
}

#endif /* SRC_NODES_SUFFIXSTORAGE_H_ */
