/*
 * NodeAddressContent.h
 *
 *  Created on: Feb 29, 2016
 *      Author: max
 */

#ifndef NODEADDRESSCONTENT_H_
#define NODEADDRESSCONTENT_H_

template <unsigned int DIM>
class Node;

template <unsigned int DIM>
struct NodeAddressContent {
	// TODO reduce bool values or align them!
	bool exists;
	bool hasSubnode;
	bool directlyStoredSuffix;
	bool hasSpecialPointer;
	int id;
	unsigned long address;
	union {
		Node<DIM>* subnode;
		const unsigned long* suffixStartBlock;
		unsigned int suffixStartBlockIndex;
		unsigned long suffix;
		uintptr_t specialPointer;
	};

	const unsigned long* getSuffixStartBlock() const {
		assert (exists && !hasSubnode);
		return (directlyStoredSuffix)? &suffix : suffixStartBlock;
	}
};


#endif /* NODEADDRESSCONTENT_H_ */
