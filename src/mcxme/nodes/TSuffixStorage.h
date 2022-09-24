/*
 * TSuffixStorage.h
 *
 *  Created on: Jun 26, 2016
 *      Author: max
 */

#ifndef SRC_NODES_TSUFFIXSTORAGE_H_
#define SRC_NODES_TSUFFIXSTORAGE_H_

template <unsigned int DIM>
class SizeVisitor;

class TSuffixStorage {
	template <unsigned int D>
	friend class SizeVisitor;
public:
	TSuffixStorage() {};
	virtual ~TSuffixStorage() {};
	virtual bool canStoreBits(size_t nBitsToStore) const =0;
	virtual unsigned int getTotalBlocksToStoreAdditionalSuffix(size_t nSuffixBits) const =0;
	// override the given index blocks with the last index blocks
	virtual unsigned int overrideBlocksWithLast(size_t nBits, unsigned int overrideStartBlockIndex) =0;
	virtual std::pair<unsigned long*, unsigned int> reserveBits(size_t nBits) =0;
	virtual void copyFrom(const TSuffixStorage& other) =0;
	virtual void clear() =0;
	virtual void clearLast(size_t nBits) =0;
	virtual unsigned int getNMaxStorageBlocks() const =0;
	virtual unsigned int getNCurrentStorageBlocks() const =0;
	virtual unsigned long getBlock(unsigned int index) const =0;
	virtual bool empty() const =0;
	virtual unsigned long* getPointerFromIndex(unsigned int index) const =0;
	virtual unsigned int getIndexFromPointer(unsigned long* pointer) const =0;
	virtual size_t getByteSize() const =0;

	size_t getNStoredSuffixes(size_t suffixBits) const;
};

size_t TSuffixStorage::getNStoredSuffixes(size_t suffixBits) const {
	const size_t currentBlocks = this->getNCurrentStorageBlocks();
	const size_t blocksPerSuffix = 1 + (suffixBits - 1) / (8 * sizeof (unsigned long));
	assert (currentBlocks % blocksPerSuffix == 0);
	return currentBlocks / blocksPerSuffix;
}

#endif /* SRC_NODES_TSUFFIXSTORAGE_H_ */
