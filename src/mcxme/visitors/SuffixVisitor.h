/*
 * SuffixVisitor.h
 *
 *  Created on: May 19, 2016
 *      Author: max
 */

#ifndef SRC_VISITORS_SUFFIXVISITOR_H_
#define SRC_VISITORS_SUFFIXVISITOR_H_


#include "../visitors/Visitor.h"

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class TNode;

template <unsigned int DIM>
class SuffixVisitor: public Visitor<DIM> {
	template <unsigned int D>
	friend std::ostream& operator <<(std::ostream &out, const SuffixVisitor<D>& v);
public:
	SuffixVisitor();
	virtual ~SuffixVisitor();

	template <unsigned int WIDTH>
	void visitSub(PHTree<DIM, WIDTH>* tree);
	template <unsigned int PREF_BLOCKS, unsigned int N>
	void visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth, unsigned int index);
	template <unsigned int PREF_BLOCKS>
	void visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth, unsigned int index);
	virtual void reset() override;

	unsigned long getPrefixSharedBits() const;
	unsigned long getPrefixBitsWithoutSharing() const;

protected:
	std::ostream& output(std::ostream &out) const override;

private:
	template <unsigned int PREF_BLOCKS>
	void visitGeneral(const TNode<DIM, PREF_BLOCKS>* node, unsigned int index);

	unsigned long externallyStoredSuffixBits;
	unsigned long externalSuffixBlocks;
	unsigned long internallyStoredSuffixBits;
	unsigned long internallyStoredSuffixes;
	unsigned long visitedSuffixes;
	unsigned long treeWidth;
	unsigned long totalUnusedLeafBlocks;
	unsigned long totalLeaves;
	vector<unsigned int> suffixBlockHistogram;
	vector<unsigned int> leafMaxSizeHistogram;
	vector<unsigned int> leafFilledSizeHistogram;
};

#include "../nodes/TNode.h"
#include "../iterators/NodeIterator.h"

template <unsigned int DIM>
SuffixVisitor<DIM>::SuffixVisitor() : Visitor<DIM>() {
	reset();
}

template <unsigned int DIM>
SuffixVisitor<DIM>::~SuffixVisitor() {
}

template <unsigned int DIM>
template <unsigned int WIDTH>
void SuffixVisitor<DIM>::visitSub(PHTree<DIM, WIDTH>* tree) {
	treeWidth = WIDTH;
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS, unsigned int N>
void SuffixVisitor<DIM>::visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int , unsigned int index) {
	this->template visitGeneral<PREF_BLOCKS>(node, index);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void SuffixVisitor<DIM>::visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int , unsigned int index) {
	this->template visitGeneral<PREF_BLOCKS>(node, index);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void SuffixVisitor<DIM>::visitGeneral(const TNode<DIM, PREF_BLOCKS>* node, unsigned int index) {

	const TSuffixStorage* storage = node->getSuffixStorage();
	const unsigned int leafSize = (storage)? storage->getNMaxStorageBlocks() : 0;
	const unsigned int filledSize = (storage)? storage->getNCurrentStorageBlocks() : 0;
	assert (leafSize >= filledSize);
	totalUnusedLeafBlocks += leafSize - filledSize;
	if (leafSize > 0) {++totalLeaves;}
	while (leafSize >= leafMaxSizeHistogram.size()) {
		leafMaxSizeHistogram.push_back(0);
		leafFilledSizeHistogram.push_back(0);
	}

	++leafMaxSizeHistogram[leafSize];
	++leafFilledSizeHistogram[filledSize];

	const unsigned long currentSuffixBits = DIM * (treeWidth - 1 - index - 1);
	assert (currentSuffixBits <= DIM * (treeWidth - 1));

	NodeIterator<DIM>* it;
	NodeIterator<DIM>* endIt = node->end();
	for (it = node->begin(); (*it) != *endIt; ++(*it)) {
		NodeAddressContent<DIM> content = *(*it);
		if (!content.hasSubnode) {
			const unsigned int suffixBlocks = 1u + currentSuffixBits / MultiDimBitset<DIM>::bitsPerBlock;
			while (suffixBlocks >= suffixBlockHistogram.size()) {
				suffixBlockHistogram.push_back(0);
			}
			++suffixBlockHistogram[suffixBlocks];
			++visitedSuffixes;
			if (content.directlyStoredSuffix) {
				++internallyStoredSuffixes;
				assert (currentSuffixBits <= sizeof (uintptr_t) * 8 - 2);
				internallyStoredSuffixBits += currentSuffixBits;
			} else {
				externallyStoredSuffixBits += currentSuffixBits;
				externalSuffixBlocks += suffixBlocks;
			}
		}
	}

	delete it;
	delete endIt;
}

template <unsigned int DIM>
void SuffixVisitor<DIM>::reset() {
	visitedSuffixes = 0;
	internallyStoredSuffixes = 0;
	treeWidth = 0;
	externallyStoredSuffixBits = 0;
	externalSuffixBlocks = 0;
	internallyStoredSuffixBits = 0;
	totalUnusedLeafBlocks = 0;
	totalLeaves = 0;
	suffixBlockHistogram.clear();
	leafMaxSizeHistogram.clear();
	leafFilledSizeHistogram.clear();
}

template <unsigned int D>
std::ostream& operator <<(std::ostream &out, const SuffixVisitor<D>& v) {
	return v.output(out);
}

template <unsigned int DIM>
std::ostream& SuffixVisitor<DIM>::output(std::ostream &out) const {

	const size_t externallyStoredSuffixes = visitedSuffixes - internallyStoredSuffixes;
	const double internallyStoredRatioPercent = double(internallyStoredSuffixes) / double(visitedSuffixes) * 100.0;
	const double avgInternalSuffixBits = double(internallyStoredSuffixBits) / double(internallyStoredSuffixes);
	const double avgExternalSuffixBits = double(externallyStoredSuffixBits) / double(externallyStoredSuffixes);
	const double avgExternalSuffixBlocks = double(externalSuffixBlocks) / double(externallyStoredSuffixes);

	out << "suffixes internally stored: " << internallyStoredSuffixes << " / "
				<< internallyStoredRatioPercent << "% (avg "
				<< avgInternalSuffixBits << " bits), avg external bits: " << avgExternalSuffixBits << " bits ("
				<< avgExternalSuffixBlocks << " blocks)" << endl;
	out << "suffix block histogram:" << endl;
	for (unsigned i = 0; i < suffixBlockHistogram.size(); ++i) {
		out << "\t" << i << " block(s): " << suffixBlockHistogram[i] << endl;
	}
	out << "unused leaf storage: " << float(totalUnusedLeafBlocks) / totalLeaves << " block(s) per leaf" << endl;
	out << "leaf size histogram:" << endl;
	assert (leafMaxSizeHistogram.size() == leafFilledSizeHistogram.size());
	out << "\t\tcapacity |\tactually used" << endl;
	for (unsigned i = 0; i < leafMaxSizeHistogram.size(); ++i) {
		out << "\t" << i << " block(s): " << leafMaxSizeHistogram[i] << " |\t" << leafFilledSizeHistogram[i] << endl;
	}

	return out;
}


#endif /* SRC_VISITORS_SUFFIXVISITOR_H_ */
