/*
 * MultiDimBitset.h
 *
 *  Created on: Apr 18, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_MULTIDIMBITSET_H_
#define SRC_UTIL_MULTIDIMBITSET_H_

#include <vector>
#include <iostream>

template <unsigned int DIM>
class MultiDimBitset {
public:
	static const size_t bitsPerBlock = sizeof(unsigned long) * 8;

	template <unsigned int WIDTH>
	static void toBitset(const std::vector<unsigned long> &values, unsigned long* const outStartBlock);

	static std::pair<bool, size_t> compare(const unsigned long* const startBlock, unsigned int nBits,
			size_t fromIndex, size_t toIndex, const unsigned long* const otherStartBlock, unsigned int otherNBits);

	static size_t compareFullAligned(const unsigned long* const startBlock,
			unsigned int nBits, const unsigned long* const otherStartBlock);

	static std::vector<unsigned long> toLongs(const unsigned long* const fromStartBlock, size_t nBits);

	static unsigned long interleaveBits(const unsigned long* const fromStartBlock, size_t index, size_t nBits);

	// TODO add: remove x <= 64 bits within block method by applying a bit mask
	static void removeHighestBits(const unsigned long* const startBlock, unsigned int nBits, size_t nBitsToRemove, unsigned long* const outStartBlock);

	static void removeHighestBits(unsigned long* const startBlock, unsigned int nBits, size_t nBitsToRemove);

	static void duplicateHighestBits(const unsigned long* startBlock, unsigned int nBits, unsigned int nBitsToDuplicate, unsigned long* const outStartBlock);

	static void duplicateBits(const unsigned long* startBlock, unsigned int nBits,
			unsigned int nBitsToDuplicate, unsigned long* const outStartBlock);

	static void duplicateLowestBitsAligned(const unsigned long* startBlock,
			unsigned nBitsToCopy, unsigned long* const outStartBlock);

	static size_t calculateLongestCommonPrefix(const unsigned long* const startBlock, unsigned int nBits, size_t startIndex,
			const unsigned long* const otherStartBlock, unsigned int otherNBits, unsigned long* const outStartBlock);

	static void pushBackValue(unsigned long value, unsigned long* pushStartBlock, unsigned int nBits);

	static void clearValue(unsigned long* startBlock, unsigned int lsbBits);

	static void pushBackBitset(const unsigned long* fromStartBlock, unsigned int fromNBits,
			unsigned long* const pushToStartBlock, unsigned int toNBits);

	static std::pair<unsigned int, unsigned int> compareSmallerEqual(const unsigned long* v1Start,
			const unsigned long* v2Start, unsigned int nBits, unsigned int skipLowestNBits, unsigned int equalDims);

	static bool checkRangeUnset(const unsigned long* startBlock, unsigned int nBits,
			unsigned int lsbStartBitIndex);

	static std::ostream& output(std::ostream &os, const unsigned long* const startBlock, unsigned int nBits);
private:
	static inline std::pair<bool, size_t> compareAlignedBlocks(const unsigned long b1, const unsigned long b2);

	static const unsigned long filledBlock = -1; // a block with all bits set to 1
	// TODO reuse this mask in all the methods
	static const unsigned long addressMask = (1uL << DIM) - 1;

	//    msb         lsb
	// off|   DIM DIM DIM
	//   0[01 001 001 001]
	static unsigned long dimBlock;
	static const unsigned int dimBlockOffset = (bitsPerBlock % DIM == 0)? 0
			: DIM - (bitsPerBlock - (((bitsPerBlock / DIM) * DIM) % bitsPerBlock));
	static unsigned long initDimBlock();
};

template <unsigned int DIM>
unsigned long MultiDimBitset<DIM>::dimBlock = initDimBlock();

#include <string.h>
#include <float.h>
#include <math.h>
#include "../libmorton/include/morton.h"

using namespace std;

template <unsigned int DIM>
unsigned long MultiDimBitset<DIM>::initDimBlock() {
	assert (dimBlockOffset < DIM);

	const size_t dimChunksPerBlock = 1 + bitsPerBlock / DIM;
	unsigned long dimBlock = 0;
	for (unsigned i = 0; i < dimChunksPerBlock; ++i) {
		dimBlock <<= DIM;
		dimBlock |= 1uL;
	}

	return dimBlock;
}

template <unsigned int DIM>
template <unsigned int WIDTH>
void MultiDimBitset<DIM>::toBitset(const std::vector<unsigned long> &values, unsigned long* outStartBlock) {
	//     example 2 Dim, 8 Bit: (    10   ,     5    )
	//    binary representation: (0000 1010, 0000 0101)
	//  				  index:    12    8    4    0
	// after conversion (mixed): (0000 0000 0110 0110) = 102
	// first dimension (mask)  : (0101 0101 0101 0101)
	// second dimension (mask) : (1010 1010 1010 1010)

	assert (sizeof (unsigned long) * 8 >= WIDTH);
	assert (values.size() == DIM);
	assert (outStartBlock[0] == 0);

	if (DIM == 2) {
		// want: (a1 b1) (a2 b2) ...
		const unsigned long v1 = values[0];
		const unsigned long v2 = values[1];
		if (WIDTH <= 32) {
			(*outStartBlock) = morton2D_64_encode(v1, v2);
		} else {
			*outStartBlock = morton2D_64_encode(v1, v2);
			*(outStartBlock + 1) = morton2D_64_encode(v1 >> 32, v2 >> 32);
		}
	} else if (DIM == 3) {
		//       [       first    block       ][ second block ...
		// want: (a1 b1 c1) (a2 b2 c2) ... (ai, bi, ci)
		// need to change order in each block because 64 % 3 == 1
		const unsigned long v1 = values[0];
		const unsigned long v2 = values[1];
		const unsigned long v3 = values[2];
		const unsigned long inter1 = morton3D_64_encode(v1, v2, v3);
		const unsigned long inter2 = morton3D_64_encode(v2 >> 21, v3 >> 21, v1 >> 22);
		const unsigned long inter3 = morton3D_64_encode(v3 >> 42, v2 >> 43, v1 >> 43);
		*outStartBlock = inter1;
		if (WIDTH > 21) {
			*(outStartBlock + 1) = inter2;
		}
		if (WIDTH > 42) {
			*(outStartBlock + 2) = inter3;
		}
	} else if (DIM == 4) {
		// want: (a1 b1 c1 d1) (a2 b2 c2 d2) ...
		// need to do: morton (morton(a, c), morton(b, d))
		const unsigned long v1 = values[0];
		const unsigned long v2 = values[1];
		const unsigned long v3 = values[2];
		const unsigned long v4 = values[3];
		const unsigned long tmp1Inter1 = morton2D_64_encode(v1, v3);
		const unsigned long tmp1Inter2 = morton2D_64_encode(v1 >> 32, v3 >> 32);
		const unsigned long tmp2Inter1 = morton2D_64_encode(v2, v4);
		const unsigned long tmp2Inter2 = morton2D_64_encode(v2 >> 32, v4 >> 32);
		const unsigned long inter1 = morton2D_64_encode(tmp1Inter1, tmp2Inter1);
		const unsigned long inter2 = morton2D_64_encode(tmp1Inter1 >> 32, tmp2Inter1 >> 32);
		const unsigned long inter3 = morton2D_64_encode(tmp1Inter2, tmp2Inter2);
		const unsigned long inter4 = morton2D_64_encode(tmp1Inter2 >> 32, tmp2Inter2 >> 32);
		*(outStartBlock) = inter1;
		if (WIDTH > 16) {
			*(outStartBlock + 1) = inter2;
		}
		if (WIDTH > 32) {
			*(outStartBlock + 2) = inter3;
		}
		if (WIDTH > 48) {
			*(outStartBlock + 3) = inter4;
		}
	} else if (DIM == 6) {
		// want: (a1 b1 c1 d1 e1 f1) (a2 b2 c2 d2 e2 f2) ...
		// need to do: morton (morton(a, c, e), morton(b, d, f))
		const unsigned long v1 = values[0];
		const unsigned long v2 = values[1];
		const unsigned long v3 = values[2];
		const unsigned long v4 = values[3];
		const unsigned long v5 = values[4];
		const unsigned long v6 = values[5];
		const unsigned long tmp1Inter1 = morton3D_64_encode(v1, v3, v5);
		const unsigned long tmp1Inter2 = morton3D_64_encode(v3 >> 21, v5 >> 21, v1 >> 22);
		const unsigned long tmp1Inter3 = morton3D_64_encode(v5 >> 42, v1 >> 43, v3 >> 43);
		const unsigned long tmp2Inter1 = morton3D_64_encode(v2, v4, v6);
		const unsigned long tmp2Inter2 = morton3D_64_encode(v4 >> 21, v6 >> 21, v2 >> 22);
		const unsigned long tmp2Inter3 = morton3D_64_encode(v6 >> 42, v2 >> 43, v4 >> 43);
		const unsigned long inter1 = morton2D_64_encode(tmp1Inter1, tmp2Inter1);
		const unsigned long inter2 = morton2D_64_encode(tmp1Inter1 >> 32, tmp2Inter1 >> 32);
		const unsigned long inter3 = morton2D_64_encode(tmp1Inter2, tmp2Inter2);
		const unsigned long inter4 = morton2D_64_encode(tmp1Inter2 >> 32, tmp2Inter2 >> 32);
		const unsigned long inter5 = morton2D_64_encode(tmp1Inter3, tmp2Inter3);
		const unsigned long inter6 = morton2D_64_encode(tmp1Inter3 >> 32, tmp2Inter3 >> 32);
		*(outStartBlock) = inter1;
		if (WIDTH > 10) {
			*(outStartBlock + 1) = inter2;
		}
		if (WIDTH > 21) {
			*(outStartBlock + 2) = inter3;
		}
		if (WIDTH > 32) {
			*(outStartBlock + 3) = inter4;
		}
		if (WIDTH > 42) {
			*(outStartBlock + 4) = inter5;
		}
		if (WIDTH > 53) {
			*(outStartBlock + 5) = inter6;
		}
	} else {
		// TODO if BMI2 is available the operation PEXT should be used!
		for (size_t d = 0; d < DIM; ++d) {
			assert (WIDTH == 64 || values[d] <= (1uL << WIDTH) - 1uL);
			const unsigned long value = values[d];
			for (size_t i = 0; i < WIDTH; ++i) {
				// extract j-th least significant bit from int
				// TODO use CPU instruction for interleaving
				const unsigned long bit = (value >> i) & 1;
				const auto block = (DIM * i + d) / bitsPerBlock;
				const auto index = (DIM * i + d) % bitsPerBlock;
				*(outStartBlock + block) |= bit << index;
			}
		}
	}

//	assert(toLongs(outStartBlock, DIM * WIDTH) == values); // TODO TZ fix this?!?
}

template <unsigned int DIM>
vector<unsigned long> MultiDimBitset<DIM>::toLongs(const unsigned long* fromStartBlock, size_t nBits) {

		vector<unsigned long> numericalValues(DIM, 0);
		const size_t bitLength = nBits / DIM;
		// TODO change loop order?
		for (size_t i = 0; i < bitLength; ++i) {
			for (size_t d = 0; d < DIM; ++d) {
				const size_t blockIndex = (DIM * i + d) / bitsPerBlock;
				const size_t bitIndex = (DIM * i + d) % bitsPerBlock;
				const size_t block = fromStartBlock[blockIndex];
				const size_t bit = (block >> bitIndex) & 1;
				numericalValues[d] += bit << i;
			}
		}

		return numericalValues;
}

template <unsigned int DIM>
std::ostream& MultiDimBitset<DIM>::output(std::ostream &os, const unsigned long* startBlock, unsigned int nBits) {
	// bits are stored in ascending order so they need to be printed in reverse oder
	assert (nBits % DIM == 0);
	const size_t bitLength = nBits / DIM;
	os << "(";
	for (size_t d = 0; d < DIM; ++d) {
		os << "(";
		for (long i = bitLength - 1; i >= 0; --i) {
			const auto blockIndex = (DIM * i + d) / bitsPerBlock;
			const auto bitIndex = (DIM * i + d) % bitsPerBlock;
			const auto block = *(startBlock + blockIndex);
			const auto bit = ((block & (1uL << bitIndex)) >> bitIndex);
			os << bit;
		}
		os << ")";
		if (d < DIM - 1) {
			os << ", ";
		}
	}

	os << ")";
	return os;
}

template <unsigned int DIM>

// fromMsbIndex is the first bit in the comparison
// toMsbIndex is the first bit after the comparison
// returns the length of the longest common prefix
std::pair<bool, size_t> MultiDimBitset<DIM>::compare(const unsigned long* startBlock, unsigned int nBits,
		size_t fromMsbIndex, size_t toMsbIndex, const unsigned long* otherStartBlock, unsigned int otherNBits) {
	assert (DIM > 0 && nBits % DIM == 0);
	assert (fromMsbIndex < toMsbIndex && toMsbIndex <= nBits / DIM);
	const size_t upper = nBits - DIM * fromMsbIndex;
	const size_t lower = nBits - DIM * toMsbIndex;
	const size_t b_max = bitsPerBlock;
	assert(lower < upper && upper < nBits);
	assert(upper - lower == otherNBits);

	// aligns blocks from the given block series to the other block series and compares blockwise
	//   (lower) |          | (upper)
	// other:    [    ][    ]
	// this : [  | ][    ][ |  ]
	//        -y-   -y-   -x-
	const size_t x = upper % b_max;
	const size_t y = lower % b_max;
	const long block_high = upper / b_max;
	const long block_low = lower / b_max;
	// extract lower part of a block [  **y**  |  b_max - y  ]
	const unsigned long lowerBlockMask = (1uL << y) - 1uL;
	// extract upper part of a block [  y  |  **b_max - y**  ]
	const unsigned long upperBlockMask = filledBlock << y;
	// removes bits after upper index from highest block: [  x  | 0 ... 0]
	const unsigned long highestBlockMask = (1uL << x) - 1uL;
	assert ((lowerBlockMask ^ upperBlockMask) == filledBlock);

	long longestCommonPrefix = 0;
	bool allDimSame = true;
	unsigned long currentAlignedBlock = 0uL;
	const unsigned long highestBlock = startBlock[block_high] & highestBlockMask;
	// handle the highest block depending on the given offsets
	if (x > y) {
		// the x range overlaps with the necessary slot so split it
		// from: [        x        |     ]
		// into: [  y  |  rest     | free]
		// compare [rest | free] to highest block and remove the rest
		const unsigned long highestAlignedBlock = highestBlock >> y;
		const size_t otherHighestBlockIndex = (otherNBits - 1) / bitsPerBlock;
		const unsigned long otherHighestBlock = otherStartBlock[otherHighestBlockIndex];
		pair<bool, size_t> comp = compareAlignedBlocks(highestAlignedBlock, otherHighestBlock);
		allDimSame = comp.first;
		const size_t highestFreeBits = b_max - x;
		assert (comp.second >= y + highestFreeBits);
		longestCommonPrefix = comp.second - y - highestFreeBits;
		// [  y  |    free   ]
		currentAlignedBlock = highestBlock & lowerBlockMask;
	} else {
		// attention: this cannot happen if the range is within a single block
		// will compare the x bits of this block in the following for-loop
		currentAlignedBlock = highestBlock;
		// will later be set to: [ b_max - y  |  x  |  free]
		// set the prefix length to the number of trailing free bits
		longestCommonPrefix = -(y - x);
	}

	for (long block = block_high - 1; block >= block_low && allDimSame; --block) {
		// [  y(t)  |  b_max - y  (t)]
		unsigned long thisBlock = startBlock[block];
		// [  y(t)  |  0    ...     0]
		const unsigned long lowerThisBlockPart = thisBlock & lowerBlockMask;
		// [0 ... 0 |  b_max - y  (t)]
		const unsigned long upperThisBlockPart = thisBlock & upperBlockMask;
		// [ 0    ...    0 |  y(t+1) ]
		currentAlignedBlock <<= b_max - y;
		// [ b_max - y  (t)|  y(t+1) ]
		currentAlignedBlock |= (upperThisBlockPart >> y);

		assert (0 <= block - block_low && otherNBits > block - block_low);
		const unsigned long otherBlock = otherStartBlock[block - block_low];
		pair<bool, size_t> comp = compareAlignedBlocks(currentAlignedBlock, otherBlock);
		allDimSame = comp.first;
		longestCommonPrefix += comp.second;

		// [  y(t)  |  0    ...     0]
		currentAlignedBlock = lowerThisBlockPart;
	}

	longestCommonPrefix /= DIM;
	assert (longestCommonPrefix >= 0 && longestCommonPrefix <= long(toMsbIndex - fromMsbIndex));
	assert (!allDimSame || otherNBits / DIM == longestCommonPrefix);
	return pair<bool, size_t>(allDimSame, longestCommonPrefix);
}

template <unsigned int DIM>
size_t MultiDimBitset<DIM>::compareFullAligned(const unsigned long* const startBlock,
			unsigned int nBits, const unsigned long* const otherStartBlock) {

	const size_t blockHighIndex = (nBits - 1uL) / bitsPerBlock;
	size_t nCommonPrefixBits = 0;
	for (unsigned block = blockHighIndex; block != (-1u); --block) {
		const pair<bool, size_t> comp = compareAlignedBlocks(startBlock[block], otherStartBlock[block]);
		nCommonPrefixBits += comp.second;
		if (!comp.first) { break; }
	}

	// remove the highest bits that were not set in both highest blocks from the count
	const size_t nCommonUnsetBits = ((blockHighIndex + 1) * bitsPerBlock) - nBits;
	nCommonPrefixBits = (nCommonPrefixBits - nCommonUnsetBits) / DIM;
	return nCommonPrefixBits;
}

template <unsigned int DIM>
pair<bool, size_t> MultiDimBitset<DIM>::compareAlignedBlocks(const unsigned long b1, const unsigned long b2) {
	const unsigned long comparison = b1 ^ b2;
	size_t longestCommonPrefix = sizeof (unsigned long) * 8;
	const bool allDimSame = comparison == 0;
	if (!allDimSame) {
		longestCommonPrefix = __builtin_clzl(comparison);
	}

	return pair<bool, size_t>(allDimSame, longestCommonPrefix);
}

template <unsigned int DIM>
unsigned long MultiDimBitset<DIM>::interleaveBits(const unsigned long* fromStartBlock, size_t msbIndex, size_t nBits) {
		assert (nBits >= DIM * (msbIndex + 1));
		assert (nBits % DIM == 0);
		assert (nBits / DIM > msbIndex);

		// TODO use blockwise operation
		unsigned long hcAddress;
		// the index addresses bits in descending order but internally bits are stored ascendingly
		const size_t lsbStartIndex = nBits - DIM * (msbIndex + 1);
		const size_t lsbEndIndex = lsbStartIndex + DIM;
		const size_t b_max = bitsPerBlock;
		const size_t startBlock = lsbStartIndex / b_max;
		const size_t endBlock = lsbEndIndex / b_max;
		assert (endBlock - startBlock <= 1);

		const unsigned long block = fromStartBlock[startBlock];
		const size_t lower = lsbStartIndex % b_max;

		if (startBlock == endBlock || lsbEndIndex % b_max == 0) {
			const unsigned long maskForInterleaved = (1uL << DIM) - 1uL; // TODO can be precalculated an reused
			hcAddress = (block >> lower) & maskForInterleaved;
		} else {
			const size_t upper = lsbEndIndex % b_max;
			const unsigned long maskForStartBlock = filledBlock << lower;
			const unsigned long maskForEndBlock = (1uL << upper) - 1;
			const unsigned long secondBlock = fromStartBlock[endBlock];
			const unsigned long shifSecond = (startBlock + 1) * b_max - lsbStartIndex;
			assert (shifSecond < DIM && shifSecond > 0);
			hcAddress = ((block & maskForStartBlock) >> lower) | ((secondBlock & maskForEndBlock) << shifSecond);
		}

		assert (hcAddress < (1ul<<DIM));
		return hcAddress;
}



template <unsigned int DIM>
void MultiDimBitset<DIM>::removeHighestBits(const unsigned long* startBlock,
		unsigned int nBits, size_t nBitsToRemove, unsigned long* outStartBlock) {
	assert (outStartBlock);
	assert (nBits >= nBitsToRemove * DIM);
	assert (nBitsToRemove > 0);

	const auto initialSize = nBits;
	if (nBitsToRemove * DIM == initialSize) {
		// no need to remove all bits
		return;
	}

	// copy all blocks before the block to be cut
	const auto remainingBits = initialSize - nBitsToRemove * DIM;
	const auto cutBlockIndex = remainingBits / bitsPerBlock;
	for (size_t b = 0; b < cutBlockIndex; ++b) {
		const unsigned long block = startBlock[b];
		outStartBlock[b] = block;
	}

	// copy the block in which the cut occurred and cut highest bits
	const auto cutBits = remainingBits % bitsPerBlock;
	if (cutBits != 0) {
		const auto lowestBitsMask = filledBlock >> (bitsPerBlock - cutBits);
		const unsigned long block = startBlock[cutBlockIndex] & lowestBitsMask;
		outStartBlock[cutBlockIndex] = block;
	}

	// removed further bits by simply not copying them
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::removeHighestBits(unsigned long* startBlock,
		unsigned int nBits, size_t nBitsToRemove) {
	assert (startBlock);
	assert (nBits >= nBitsToRemove && (nBitsToRemove > 0) && (nBitsToRemove % DIM == 0));

	// lsb                        msb
	//      <-remove-><-not changed->
	// <----nBits---->
	// [part][full][part][ -- ][ -- ]

	const unsigned int firstClearedBit = nBits - nBitsToRemove;
	const unsigned int firstClearBlockIndex = firstClearedBit / bitsPerBlock;
	const unsigned int firstClearBitIndex = firstClearedBit % bitsPerBlock;
	const unsigned int lastClearBlockIndex = (nBits - 1) / bitsPerBlock;
	const unsigned int lastClearBlockBitIndex = (nBits - 1) % bitsPerBlock;

	if (firstClearBlockIndex == lastClearBlockIndex) {
		// clear within one block
		assert (nBitsToRemove == lastClearBlockBitIndex - firstClearBitIndex + 1);
		const unsigned long singleBlockMask = (nBitsToRemove == bitsPerBlock) ?
						0 :
						~(((1uL << nBitsToRemove) - 1uL) << firstClearBitIndex);
		assert (nBitsToRemove == bitsPerBlock || singleBlockMask != 0);
		startBlock[firstClearBlockIndex] &= singleBlockMask;
	} else {
		// clear within several block (lsb and msb masks can be applied!)
		// retain the lsb bits that are not within the removal scope
		if (firstClearBitIndex == 0) startBlock[firstClearBlockIndex] = 0;
		else {
			const unsigned long lsbBlockMask = filledBlock >> (bitsPerBlock - firstClearBitIndex);
			startBlock[firstClearBlockIndex] &= lsbBlockMask;
		}

		if (lastClearBlockIndex - firstClearBlockIndex > 1) {
			// fully clear all blocks that are completely in the removal range:
			const unsigned int nFullClearBlocks = lastClearBlockIndex - firstClearBlockIndex - 1;
			for (unsigned i = 0; i < nFullClearBlocks; ++i) {
				startBlock[firstClearBlockIndex + 1 + i] = 0;
			}
		}

		// retain the msb bits that are not within the removal scope
		if (lastClearBlockBitIndex == bitsPerBlock - 1) startBlock[lastClearBlockIndex] = 0;
		else {
			const unsigned long msbBlockMask = filledBlock << (lastClearBlockBitIndex + 1);
			startBlock[lastClearBlockIndex] &= msbBlockMask;
		}
	}
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::duplicateHighestBits(const unsigned long* startBlock,
		unsigned int nBits, unsigned int nBitsToDuplicate, unsigned long* outStartBlock) {
	assert (outStartBlock);
	assert (nBits >= DIM * nBitsToDuplicate);
	assert (nBitsToDuplicate > 0);

	// duplicate:       <--------------nBitsToDuplicate------------------>
	// this:      | v_1 \   r_1     | , | v_2 \   r_2     | , | v_3 \ r_3 \ free|
	//            | vL  \    rL     | , |-----\-----------| , |-----\-----------|
	// result:    |   r_1    \ v_2  | , |   r_2     \ v_3 | , |  r_3 \    free  |
	const size_t neglectedLSBBits = nBits - DIM * nBitsToDuplicate;
	const size_t startCopyBlock = neglectedLSBBits / bitsPerBlock;
	const size_t lastBlockIndex = (nBits - 1) / bitsPerBlock;
	const unsigned long vL = neglectedLSBBits % bitsPerBlock;
	const unsigned long rL = bitsPerBlock - vL;
	const unsigned long vMask = (1uL << vL) - 1;
	const unsigned long rMask = filledBlock << vL;
	assert ((vL + rL == bitsPerBlock) && ((vMask ^ rMask) == filledBlock));

	const unsigned long sourceBlock = startBlock[startCopyBlock];
	const unsigned long r = (sourceBlock & rMask) >> vL;
	unsigned long currentBlock = r;
	// add current value to free space in current block and set remainder to next block
	for (size_t b = startCopyBlock + 1; b < lastBlockIndex + 1; ++b) {
		const unsigned long sourceBlock = startBlock[b];
		const unsigned long v = (sourceBlock & vMask) << rL;
		const unsigned long r = (sourceBlock & rMask) >> vL;
		currentBlock |= v;
		assert (outStartBlock[b - (startCopyBlock + 1)] == 0);
		outStartBlock[b - (startCopyBlock + 1)] = currentBlock;
		currentBlock = r;
	}

	const size_t lastBit = nBits % bitsPerBlock;
	if (vL < lastBit || lastBit == 0) {
		// needed if:
		// - the last part (r) was not yet fully copied so the last bit(s) are not yet transfered
		// - the last bit index (!) is actually in the previous block
		//   (i.e. the last part is aligned with the end of a block)
		assert (outStartBlock[lastBlockIndex - startCopyBlock] == 0);
		outStartBlock[lastBlockIndex - startCopyBlock] = currentBlock;
	}
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::duplicateBits(const unsigned long* startBlock, unsigned int nBits,
		unsigned int nBitsToDuplicate, unsigned long* const outStartBlock) {

	// duplicate the bits as if they are the highest bits
	duplicateHighestBits(startBlock, nBits, nBitsToDuplicate, outStartBlock);

	assert (nBits >= DIM * nBitsToDuplicate);
	assert (nBitsToDuplicate > 0);

	// remove bits from the highest block that might be copied wrongly
	const size_t newNBits = DIM * nBitsToDuplicate;
	const size_t lastCopyBitIndex = (newNBits - 1) % bitsPerBlock;
	const size_t lastBlockIndex = (newNBits - 1) / bitsPerBlock;
	const unsigned long highestBlockMask = ~((-1uL) << (lastCopyBitIndex + 1));
	outStartBlock[lastBlockIndex] &= highestBlockMask;
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::duplicateLowestBitsAligned(const unsigned long* startBlock,
		unsigned nBitsToCopy, unsigned long* const outStartBlock) {
	assert (nBitsToCopy % DIM == 0);

	const size_t nFullCopyBlocks = nBitsToCopy / bitsPerBlock;
	for (unsigned b = 0; b < nFullCopyBlocks; ++b) {
		outStartBlock[b] = startBlock[b];
	}

	const size_t remainingBits = nBitsToCopy % bitsPerBlock;
	if (remainingBits > 0) {
		outStartBlock[nFullCopyBlocks] = startBlock[nFullCopyBlocks];
		// remove bits that were copied to much
		const unsigned long upperEndMask = ~(filledBlock << remainingBits);
		outStartBlock[nFullCopyBlocks] &= upperEndMask;
	}
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::clearValue(unsigned long* startBlock, unsigned int lsbBits) {
	assert (lsbBits % DIM == 0);

	const unsigned int startBlockIndex = lsbBits / bitsPerBlock;
	const unsigned int startBitIndex = lsbBits % bitsPerBlock;
	const unsigned int endBlockIndex = (lsbBits + DIM - 1) / bitsPerBlock;

	startBlock[startBlockIndex] &= ~(addressMask << startBitIndex);

	assert (endBlockIndex - startBlockIndex <= 1);
	if (startBlockIndex != endBlockIndex) {
		assert (startBitIndex != 0 && "otherwise the shift does not work");
		const size_t remainingCancelBits = DIM - (bitsPerBlock - startBitIndex);
		startBlock[endBlockIndex] &= filledBlock << remainingCancelBits;
	}
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::pushBackValue(unsigned long interleavedValue, unsigned long* startBlock, unsigned int nBits) {
	assert (startBlock && nBits % DIM == 0 && interleavedValue < 1<<DIM);

	// the first bit of the interleaved address ought to be placed at index nBits
	const unsigned int startBlockIndex = nBits / bitsPerBlock;
	const unsigned int startBitIndex = nBits % bitsPerBlock;
	const unsigned int endBlockIndex = (nBits + DIM - 1u) / bitsPerBlock;

	startBlock[startBlockIndex] |= interleavedValue << startBitIndex;

	assert (endBlockIndex - startBlockIndex <= 1);
	if (startBlockIndex != endBlockIndex) {
		startBlock[endBlockIndex] |= interleavedValue >> (bitsPerBlock - startBitIndex);
	}
}

template <unsigned int DIM>
void MultiDimBitset<DIM>::pushBackBitset(const unsigned long* fromStartBlock, unsigned int fromNBits,
			unsigned long* const pushToStartBlock, unsigned int toNBits) {
	assert (fromStartBlock && pushToStartBlock && fromNBits % DIM == 0 && toNBits % DIM == 0);
	assert (fromNBits > 0);

	//       <--fromNBits-->
	// from: [    ][    ][  abc]
	//       <--toNBits--->
	// to:   [    ][    ][ 000][0000][0000][xyz][....]
	//       <--toNBits---><--fromNBits--><-unchanged->
	// to:   [ -- ][ -- ][    ][    ][    ][xyz][....]

	const unsigned int fromNBitsIndex = fromNBits - 1u;
	const unsigned int fromBitIndex = fromNBitsIndex % bitsPerBlock;
	const unsigned int startPushToBlockIndex = toNBits / bitsPerBlock;
	const unsigned int startPushToBitIndex = toNBits % bitsPerBlock;
	const unsigned int pushNBlocks = 1u + fromNBitsIndex / bitsPerBlock;
	const unsigned long firstToMask = filledBlock >> (bitsPerBlock - startPushToBitIndex);
	const unsigned long secondToMask = filledBlock ^ firstToMask;

	assert (pushNBlocks > 0);
	// handle all from blocks that are fully copied
	for (unsigned i = 0; i < pushNBlocks - 1u; ++i) {
		const unsigned long fromBlock = fromStartBlock[i];
		const unsigned long firstFromPart = fromBlock & firstToMask;
		const unsigned long secondFromPart = fromBlock & secondToMask;
		pushToStartBlock[startPushToBlockIndex + i] |= firstFromPart << startPushToBitIndex;

		if (secondToMask != 0)
			pushToStartBlock[startPushToBlockIndex + i + 1] |= secondFromPart >> startPushToBitIndex;
	}

	// handle last from block if it is partial
		// the last block should only have 0s after the last known bit
		// therefore a the entire block can be shifted
		const unsigned long lastBlockFromMask = filledBlock >> (bitsPerBlock - fromBitIndex - 1);
		const unsigned long fromBlock = fromStartBlock[pushNBlocks - 1] & lastBlockFromMask;
		pushToStartBlock[startPushToBlockIndex + pushNBlocks - 1] |= fromBlock << startPushToBitIndex;
		if (secondToMask != 0 && fromBitIndex + startPushToBitIndex >= bitsPerBlock) {
			const unsigned int remainingPushBits = bitsPerBlock - startPushToBitIndex;
			assert (remainingPushBits > 0);
			pushToStartBlock[startPushToBlockIndex + pushNBlocks] |= fromBlock >> remainingPushBits;
		}

}

template <unsigned int DIM>
bool MultiDimBitset<DIM>::checkRangeUnset(const unsigned long* startBlock, unsigned int nBits,
			unsigned int lsbStartBitIndex) {
	if (nBits == 0 || nBits == lsbStartBitIndex) return true;

	assert (lsbStartBitIndex < nBits);
	assert (lsbStartBitIndex % DIM == 0);

	const unsigned int lsbEndBitIndex = nBits - 1;
	const unsigned int startCompareBlockIndex = lsbStartBitIndex / bitsPerBlock;
	const unsigned int startCompareBitIndex = lsbStartBitIndex % bitsPerBlock;
	const unsigned int endCompareBlockIndex = lsbEndBitIndex / bitsPerBlock;
	const unsigned int endCompareBitIndex = lsbEndBitIndex % bitsPerBlock;

	bool rangeBitsUnset = true;
	for (unsigned blockIndex = startCompareBlockIndex;
			blockIndex < endCompareBlockIndex + 1uL && rangeBitsUnset;
			++blockIndex) {
		unsigned long block = startBlock[blockIndex];
		if (blockIndex == startCompareBlockIndex) {
			// mask out lower bits
			block &= filledBlock << startCompareBitIndex;
		}
		if (blockIndex == endCompareBlockIndex) {
			// mask out upper bits
			block &= filledBlock >> (bitsPerBlock - endCompareBitIndex - 1u);
		}

		rangeBitsUnset = (block == 0uL);
	}

	return rangeBitsUnset;
}

template <unsigned int DIM>
pair<unsigned int, unsigned int> MultiDimBitset<DIM>::compareSmallerEqual(const unsigned long* v1Start,
		const unsigned long* v2Start, unsigned int nBits, unsigned int skipLowestNBits, unsigned int equalDims) {
	assert (nBits > skipLowestNBits && nBits > 0);
	assert (nBits % DIM == 0 && skipLowestNBits % DIM == 0);
	assert (equalDims > 0 && "otherwise there is no reason to call this method");

	bool isSmaller[DIM] = {};
	bool isEqual[DIM] = {};
	bool allCompared = false;
	for (unsigned d = 0; d < DIM; ++d) {
		isEqual[d] = (equalDims >> d) & 1;
	}

	const unsigned int nBitsIndex = nBits - 1u;
	const unsigned int highestBlock = nBitsIndex / bitsPerBlock;
	const unsigned int lowestBlock = skipLowestNBits / bitsPerBlock;
	assert (highestBlock >= lowestBlock);
	const unsigned int cancelLowBitIndex = skipLowestNBits % bitsPerBlock;
	const unsigned long cancelLowestBitsMask = filledBlock << cancelLowBitIndex;
	const unsigned int cancelHighBitIndex = nBitsIndex % bitsPerBlock;
	const unsigned long cancelHighestBitsMask = filledBlock >> (bitsPerBlock - cancelHighBitIndex - 1);
	assert (cancelHighestBitsMask != 0);

	for (unsigned int block = highestBlock; !allCompared; --block) {
		// apply offset to lowest dim bit mask
		const unsigned int blockOffset = (dimBlockOffset * block) % DIM;
		assert (blockOffset < DIM);

		allCompared = true;
		const unsigned long v1Block = v1Start[block];
		const unsigned long v2Block = v2Start[block];
		for (unsigned d = 0; d < DIM; ++d) {
			// (in lower blocks) the result can only be overwritten if the values were equal so far
			if (!isEqual[d]) continue;

			const unsigned int dimOffset = (d + blockOffset) % DIM;
			unsigned long dimExtractionMask = dimBlock << dimOffset;
			if (block == lowestBlock) { dimExtractionMask &= cancelLowestBitsMask; }
			if (block == highestBlock) { dimExtractionMask &= cancelHighestBitsMask; }

			const unsigned long v1BlockDimExtracted = v1Block & dimExtractionMask;
			const unsigned long v2BlockDimExtracted = v2Block & dimExtractionMask;
			isSmaller[d] = v1BlockDimExtracted < v2BlockDimExtracted;
			isEqual[d] = !isSmaller[d] && (v1BlockDimExtracted == v2BlockDimExtracted);
			allCompared &= !isEqual[d];
		}

		if (block == lowestBlock) break;
	}


	unsigned int dimLowerComparison = 0;
	unsigned int dimEqualComparison = 0;
	for (unsigned d = 0; d < DIM; ++d) {
		assert (!(isSmaller[d] && isEqual[d]));
		if (isSmaller[d]) {
			dimLowerComparison |= 1uL << d;
		} else if (isEqual[d]) {
			dimEqualComparison |= 1uL << d;
		}
	}

	assert (dimLowerComparison < 1uL << DIM);
	assert (dimEqualComparison < 1uL << DIM);
	assert ((dimLowerComparison & dimEqualComparison) == 0);
	return pair<unsigned int, unsigned int>(dimLowerComparison, dimEqualComparison);
}

// find the longest common prefix starting at the msbStartIndex and sets it to the result to reference
template <unsigned int DIM>
size_t MultiDimBitset<DIM>::calculateLongestCommonPrefix(const unsigned long* startBlock,
		unsigned int nBits, size_t msbStartIndex, const unsigned long* otherStartBlock,
		unsigned int otherNBits, unsigned long* outStartBlock) {

	assert (startBlock && otherStartBlock && outStartBlock);
	assert (nBits - DIM * msbStartIndex == otherNBits);

	const size_t otherBitLength = otherNBits / DIM;
	const size_t msbEndIndex = msbStartIndex + otherBitLength;
	pair<bool, unsigned long> comparison = compare(startBlock, nBits, msbStartIndex,
			msbEndIndex, otherStartBlock, otherNBits);
	const size_t prefixLength = comparison.second;
	assert (prefixLength <= otherBitLength);
	if (prefixLength > 0) {
		duplicateHighestBits(otherStartBlock, otherNBits, prefixLength, outStartBlock);
	}

	return prefixLength;
}

#endif /* SRC_UTIL_MULTIDIMBITSET_H_ */
