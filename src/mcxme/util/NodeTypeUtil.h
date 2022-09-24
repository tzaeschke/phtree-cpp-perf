/*
 * NodeTypeUtil.h
 *
 *  Created on: Apr 7, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_NODETYPEUTIL_H_
#define SRC_UTIL_NODETYPEUTIL_H_

#include <cstdint>
#include "../nodes/LHC.h"
#include "../nodes/AHC.h"
#include "../nodes/SuffixStorage.h"
#include "../util/TEntryBuffer.h"

template <unsigned int DIM>
class Node;

template <unsigned int DIM>
class NodeTypeUtil {
public:
	template <unsigned int WIDTH>
	static Node<DIM>* buildNodeWithSuffixes(size_t prefixBits, size_t nDirectInserts, size_t nSuffixes, unsigned int suffixBits) {
		assert (suffixBits < DIM * WIDTH && (suffixBits % DIM == 0));
		assert (nSuffixes <= nDirectInserts);
		Node<DIM>* node = buildNode(prefixBits, nDirectInserts);
		if (nSuffixes > 0 && suffixBits > 0 && !node->canStoreSuffixInternally(suffixBits)) {
			const unsigned int suffixBlocks = 1 + (suffixBits - 1) / (8 * sizeof (unsigned long));
			TSuffixStorage* storage = createSuffixStorage<WIDTH>(nSuffixes * suffixBlocks);
			node->setSuffixStorage(storage);
		}

		return node;
	}

	template <unsigned int WIDTH>
	static void enlargeSuffixStorage(unsigned int suffixBlocks, Node<DIM>* node) {
		assert (suffixBlocks > 0);

		TSuffixStorage* suffixes = createSuffixStorage<WIDTH>(suffixBlocks);

		// copy the old contents
		const TSuffixStorage* oldStorage = node->getSuffixStorage();
		assert (!oldStorage || suffixes->getNMaxStorageBlocks() > oldStorage->getNMaxStorageBlocks());
		if (oldStorage) {
			suffixes->copyFrom(*oldStorage);
			delete oldStorage;
		}

		node->setSuffixStorage(suffixes);
	}

	template <unsigned int WIDTH>
	static void shrinkSuffixStorageIfPossible(Node<DIM>* node) {
		assert (node);
		const TSuffixStorage* oldStorage = node->getSuffixStorage();
		if (oldStorage) {
			const unsigned int filledSuffixBlocks = oldStorage->getNCurrentStorageBlocks();
			const unsigned int maxPosibleSuffixBlocks = oldStorage->getNMaxStorageBlocks();
			bool empty;
			const bool shrink = canShrinkSuffixStorage<WIDTH>(filledSuffixBlocks, maxPosibleSuffixBlocks, &empty);
			if (shrink) {
				TSuffixStorage* shrinkedStorage = NULL;
				if (!empty) {
					shrinkedStorage = createSuffixStorage<WIDTH>(filledSuffixBlocks);
					shrinkedStorage->copyFrom(*oldStorage);
					assert (shrinkedStorage->getNMaxStorageBlocks() < oldStorage->getNMaxStorageBlocks());
				}

				node->setSuffixStorage(shrinkedStorage);
				delete oldStorage;
			}
		}
	}

	static Node<DIM>* copyWithoutPrefix(size_t newPrefixBits, const Node<DIM>* nodeToCopy) {
		// TODO no need to duplicate suffix storage
		size_t nDirectInsert = nodeToCopy->getNumberOfContents();
		Node<DIM>* copy = buildNode(newPrefixBits, nDirectInsert);
		copyContents(*nodeToCopy, *copy);
		return copy;
	}

	static Node<DIM>* copyIntoLargerNode(size_t newNContents, const Node<DIM>* nodeToCopy) {
		// TODO make more efficient by not using iterators and a bulk insert
		const size_t prefixLength = nodeToCopy->getPrefixLength();
		Node<DIM>* copy = buildNode(prefixLength * DIM, newNContents);
		if (prefixLength > 0) {
			MultiDimBitset<DIM>::duplicateHighestBits(nodeToCopy->getFixPrefixStartBlock(),
					prefixLength * DIM, prefixLength, copy->getPrefixStartBlock());
		}

		copyContents(*nodeToCopy, *copy);
		return copy;
	}

private:

	template <unsigned int WIDTH>
	inline static bool canShrinkSuffixStorage(unsigned int newRequiredSuffixBlocks, unsigned int oldSuffixBlocks, bool* empty) {
		assert (newRequiredSuffixBlocks < oldSuffixBlocks);
		unsigned int newGrantedSuffixBlocks = 0;
		const unsigned int maxSuffixBits = (WIDTH - 1) * DIM;
		const unsigned int maxBlocksPerSuffix = 1 + (maxSuffixBits - 1) / (8 * sizeof (unsigned long));
		const unsigned int maxSuffixBlocks = maxBlocksPerSuffix * (1u << DIM);
		const float suffixRatio = float(newRequiredSuffixBlocks) / float(maxSuffixBlocks);
		if (newRequiredSuffixBlocks < 6) {
			switch (newRequiredSuffixBlocks) {
			case 0: newGrantedSuffixBlocks = 0; break;
			case 1: newGrantedSuffixBlocks = 1; break;
			case 2: newGrantedSuffixBlocks = 2; break;
			case 3: newGrantedSuffixBlocks = 3; break;
			case 4: newGrantedSuffixBlocks = 4; break;
			case 5: newGrantedSuffixBlocks = 5; break;
			default: throw runtime_error("Only supports up to 5 fixed suffix blocks right now.");
			}
		} else if (newRequiredSuffixBlocks < 8 ) {
			newGrantedSuffixBlocks = 7;
		} else if (newRequiredSuffixBlocks < 11) {
			newGrantedSuffixBlocks = 10;
		} else if (suffixRatio < 0.001) {
			newGrantedSuffixBlocks = 1 + maxSuffixBlocks / 1000;
		} else if (suffixRatio < 0.005) {
			newGrantedSuffixBlocks = 1 + 5 * maxSuffixBlocks / 1000;
		} else if (suffixRatio < 0.01) {
			newGrantedSuffixBlocks = 1 + maxSuffixBlocks / 100;
		} else if (suffixRatio < 0.05) {
			newGrantedSuffixBlocks = 1 + 5 * maxSuffixBlocks / 100;
		} else if (suffixRatio < 0.1) {
			newGrantedSuffixBlocks = 1 + 10 * maxSuffixBlocks / 100;
		} else if (suffixRatio < 0.25) {
			newGrantedSuffixBlocks = 1 + 25 * maxSuffixBlocks / 100;
		} else if (suffixRatio < 0.5) {
			newGrantedSuffixBlocks = 1 + 50 * maxSuffixBlocks / 100;
		} else if (suffixRatio < 0.75) {
			newGrantedSuffixBlocks = 1 + 75 * maxSuffixBlocks / 100;
		} else {
			newGrantedSuffixBlocks = maxSuffixBlocks;
		}

		assert (newGrantedSuffixBlocks <= oldSuffixBlocks);
		(*empty) = newGrantedSuffixBlocks == 0;
		return newGrantedSuffixBlocks != oldSuffixBlocks;
	}

	template <unsigned int WIDTH>
	inline static TSuffixStorage* createSuffixStorage(unsigned int suffixBlocks) {
		assert (suffixBlocks > 0);
		// a node has at most 2^DIM suffixes and one suffix is at most (WIDTH-1) bits long in each dimension
		const unsigned int maxSuffixBits = (WIDTH - 1) * DIM;
		const unsigned int maxBlocksPerSuffix = 1 + (maxSuffixBits - 1) / (8 * sizeof (unsigned long));
		const unsigned int maxSuffixBlocks = maxBlocksPerSuffix * (1u << DIM);
		assert (suffixBlocks <= maxSuffixBlocks);
		const float suffixRatio = float(suffixBlocks) / float(maxSuffixBlocks);
		assert (suffixRatio <= 1.0);
		TSuffixStorage* suffixes;
		if (suffixBlocks < 6) {
			switch (suffixBlocks) {
			case 1: suffixes = new SuffixStorage<1>(); break;
			case 2: suffixes = new SuffixStorage<2>(); break;
			case 3: suffixes = new SuffixStorage<3>(); break;
			case 4: suffixes = new SuffixStorage<4>(); break;
			case 5: suffixes = new SuffixStorage<5>(); break;
			default: throw runtime_error("Only supports up to 5 fixed suffix blocks right now.");
			}
		} else if (suffixBlocks < 8 ) {
			suffixes = new SuffixStorage<7>();
		} else if (suffixBlocks < 11) {
			suffixes = new SuffixStorage<10>();
		} else if (suffixRatio < 0.001) {
			suffixes = new SuffixStorage<1 + maxSuffixBlocks / 1000>();
		} else if (suffixRatio < 0.005) {
			suffixes = new SuffixStorage<1 + 5 * maxSuffixBlocks / 1000>();
		} else if (suffixRatio < 0.01) {
			suffixes = new SuffixStorage<1 + maxSuffixBlocks / 100>();
		} else if (suffixRatio < 0.05) {
			suffixes = new SuffixStorage<1 + 5 * maxSuffixBlocks / 100>();
		} else if (suffixRatio < 0.1) {
			suffixes = new SuffixStorage<1 + 10 * maxSuffixBlocks / 100>();
		} else if (suffixRatio < 0.25) {
			suffixes = new SuffixStorage<1 + 25 * maxSuffixBlocks / 100>();
		} else if (suffixRatio < 0.5) {
			suffixes = new SuffixStorage<1 + 50 * maxSuffixBlocks / 100>();
		} else if (suffixRatio < 0.75) {
			suffixes = new SuffixStorage<1 + 75 * maxSuffixBlocks / 100>();
		} else {
			suffixes = new SuffixStorage<maxSuffixBlocks>();
		}

		return suffixes;
	}

	inline static void copyContents(const Node<DIM>& from, Node<DIM>& to) {
		// copy suffixes
		assert (!to.getSuffixStorage());
		to.copySuffixStorageFrom(from);

		// TODO make more efficient by not using iterators and a bulk insert
		// copy node contents
		NodeIterator<DIM>* it = from.begin();
		it->disableResolvingSuffixIndex();
		NodeIterator<DIM>* endIt = from.end();
		for (; (*it) != *endIt; ++(*it)) {
			NodeAddressContent<DIM> content = *(*it);
			if (content.hasSubnode) {
				to.insertAtAddress(content.address, content.subnode);
			} else if (content.hasSpecialPointer) {
				to.insertAtAddress(content.address, content.specialPointer);
				TEntryBuffer<DIM>* buffer = reinterpret_cast<TEntryBuffer<DIM>*>(content.specialPointer);
				buffer->updateNode(&to);
			} else if (content.directlyStoredSuffix) {
				to.insertAtAddress(content.address, content.suffix, content.id);
			} else {
				to.insertAtAddress(content.address, content.suffixStartBlockIndex, content.id);
			}
		}

		delete it;
		delete endIt;
	}

	template <unsigned int PREF_BLOCKS>
	inline static Node<DIM>* determineNodeType(size_t prefixBits, size_t nDirectInserts) {
		assert (nDirectInserts > 0);
		const size_t prefixLength = prefixBits / DIM;
		// TODO use threshold depending on which node is smaller
		const double switchTypeAtLoadRatio = 0.75;
		if (float(nDirectInserts) / (1uL << DIM) < switchTypeAtLoadRatio) {
			return determineLhcSize<PREF_BLOCKS>(prefixLength, nDirectInserts);
		} else {
			return new AHC<DIM, PREF_BLOCKS>(prefixLength);
		}
	}

	inline static Node<DIM>* buildNode(size_t prefixBits, size_t nDirectInserts) {
			const size_t prefixBlocks = (prefixBits > 0)? 1 + ((prefixBits - 1) / (8 * sizeof (unsigned long))) : 0;
			switch (prefixBlocks) {
			case 0: return determineNodeType<0>(prefixBits, nDirectInserts);
			case 1: return determineNodeType<1>(prefixBits, nDirectInserts);
			case 2: return determineNodeType<2>(prefixBits, nDirectInserts);
			case 3: return determineNodeType<3>(prefixBits, nDirectInserts);
			case 4: return determineNodeType<4>(prefixBits, nDirectInserts);
			case 5: return determineNodeType<5>(prefixBits, nDirectInserts);
			case 6: return determineNodeType<6>(prefixBits, nDirectInserts);
			case 7: return determineNodeType<7>(prefixBits, nDirectInserts);
			case 8: return determineNodeType<8>(prefixBits, nDirectInserts);
			case 9: return determineNodeType<9>(prefixBits, nDirectInserts);
			case 10: return determineNodeType<10>(prefixBits, nDirectInserts);
			default: throw runtime_error("Only supports up to 10 prefix blocks right now.");
			}
		}

	template<unsigned int PREF_BLOCKS>
	inline static Node<DIM>* determineLhcSize(size_t prefixLength,
			size_t nDirectInserts) {

		assert(nDirectInserts > 0);
		const float insertToRatio = float(nDirectInserts) / (1u << DIM);
		assert(0 < insertToRatio && insertToRatio < 1);

		if (nDirectInserts < 3) {
			return new LHC<DIM, PREF_BLOCKS, 2>(prefixLength);
		} else if (insertToRatio < 0.1) {
			return new LHC<DIM, PREF_BLOCKS, 1 + 10 * (1 << DIM) / 100>(prefixLength);
		} else if (insertToRatio < 0.2) {
			return new LHC<DIM, PREF_BLOCKS, 1 + 20 * (1 << DIM) / 100>(prefixLength);
		} else if (insertToRatio < 0.35) {
			return new LHC<DIM, PREF_BLOCKS, 1 + 35 * (1 << DIM) / 100>(prefixLength);
		} else if (insertToRatio < 0.5) {
			return new LHC<DIM, PREF_BLOCKS, 1 + 50 * (1 << DIM) / 100>(prefixLength);
		} else if (insertToRatio < 0.75) {
			return new LHC<DIM, PREF_BLOCKS, 1 + 75 * (1 << DIM) / 100>(prefixLength);
		} else {
			return new LHC<DIM, PREF_BLOCKS, (1 << DIM)>(prefixLength);
		}
	}
};

#endif /* SRC_UTIL_NODETYPEUTIL_H_ */
