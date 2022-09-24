/*
 * PrefixSharingVisitor.h
 *
 *  Created on: Apr 25, 2016
 *      Author: max
 */

#ifndef SRC_VISITORS_PREFIXSHARINGVISITOR_H_
#define SRC_VISITORS_PREFIXSHARINGVISITOR_H_

#include "../visitors/Visitor.h"

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class TNode;

template <unsigned int DIM>
class PrefixSharingVisitor: public Visitor<DIM> {
	template <unsigned int D>
	friend std::ostream& operator <<(std::ostream &out, const PrefixSharingVisitor<D>& v);
public:
	PrefixSharingVisitor();
	virtual ~PrefixSharingVisitor();

	template <unsigned int WIDTH>
	void visitSub(PHTree<DIM, WIDTH>* tree);
	template <unsigned int PREF_BLOCKS, unsigned int N>
	void visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth);
	template <unsigned int PREF_BLOCKS>
	void visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth);
	virtual void reset() override;

	unsigned long getPrefixSharedBits() const;
	unsigned long getPrefixBitsWithoutSharing() const;

protected:
	std::ostream& output(std::ostream &out) const override;

private:
	unsigned long prefixSharedBits;
	unsigned long prefixBitsWithoutSharing;
	unsigned long nSharedPrefixes;
	vector<unsigned int> prefixBlockHistogram;

	template <unsigned int PREF_BLOCKS>
	void visitGeneral(const TNode<DIM, PREF_BLOCKS>* node);
};

#include "../nodes/TNode.h"

template <unsigned int DIM>
PrefixSharingVisitor<DIM>::PrefixSharingVisitor() : Visitor<DIM>() {
	reset();
}

template <unsigned int DIM>
PrefixSharingVisitor<DIM>::~PrefixSharingVisitor() {
}

template <unsigned int DIM>
template <unsigned int WIDTH>
void PrefixSharingVisitor<DIM>::visitSub(PHTree<DIM, WIDTH>* tree) {

}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS, unsigned int N>
void PrefixSharingVisitor<DIM>::visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth) {
	this->template visitGeneral<PREF_BLOCKS>(node);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void PrefixSharingVisitor<DIM>::visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth) {
	this->template visitGeneral<PREF_BLOCKS>(node);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void PrefixSharingVisitor<DIM>::visitGeneral(const TNode<DIM, PREF_BLOCKS>* node) {
	while (PREF_BLOCKS >= prefixBlockHistogram.size()) {
		prefixBlockHistogram.push_back(0);
	}
	++prefixBlockHistogram[PREF_BLOCKS];

	prefixSharedBits += node->getPrefixLength() * DIM;
	prefixBitsWithoutSharing += node->getPrefixLength() * DIM * node->getNumberOfContents();
	if (node->getPrefixLength() > 0) {
		++nSharedPrefixes;
	}
}

template <unsigned int DIM>
void PrefixSharingVisitor<DIM>::reset() {
	prefixBitsWithoutSharing = 0;
	prefixSharedBits = 0;
	nSharedPrefixes = 0;
	prefixBlockHistogram.clear();
}

template <unsigned int DIM>
unsigned long PrefixSharingVisitor<DIM>::getPrefixSharedBits() const {
	return prefixSharedBits;
}

template <unsigned int DIM>
unsigned long PrefixSharingVisitor<DIM>::getPrefixBitsWithoutSharing() const {
	return prefixBitsWithoutSharing;
}

template <unsigned int D>
std::ostream& operator <<(std::ostream &out, const PrefixSharingVisitor<D>& v) {
	return v.output(out);
}

template <unsigned int DIM>
std::ostream& PrefixSharingVisitor<DIM>::output(std::ostream &out) const {
	float sharingSavedPercent = float(getPrefixSharedBits()) / float(getPrefixBitsWithoutSharing());
	out << "#nodes with prefix=" << nSharedPrefixes << ", total shared prefix bits: " << getPrefixSharedBits() << " bits"
			<< " (only " << sharingSavedPercent << "% of otherwise "
			<< getPrefixBitsWithoutSharing() << " bits)" << endl;
	out << "prefix blocks histogram:" << endl;
	for (unsigned i = 0; i < prefixBlockHistogram.size(); ++i) {
		out << "\t" << i << " block(s): " << prefixBlockHistogram[i] << endl;
	}

	return out;
}

#endif /* SRC_VISITORS_PREFIXSHARINGVISITOR_H_ */
