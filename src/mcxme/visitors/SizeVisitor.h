/*
 * SizeVisitor.h
 *
 *  Created on: Apr 13, 2016
 *      Author: max
 */

#ifndef SRC_VISITORS_SIZEVISITOR_H_
#define SRC_VISITORS_SIZEVISITOR_H_

#include "../visitors/Visitor.h"
#include "../util/MultiDimBitset.h"

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class TNode;

template <unsigned int DIM>
class SizeVisitor: public Visitor<DIM> {
	template <unsigned int D>
	friend std::ostream& operator <<(std::ostream &out, const SizeVisitor<D>& v);
public:
	SizeVisitor();
	virtual ~SizeVisitor();

	template <unsigned int WIDTH>
	void visitSub(PHTree<DIM, WIDTH>* tree);
	template <unsigned int PREF_BLOCKS, unsigned int N>
	void visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth);
	template <unsigned int PREF_BLOCKS>
	void visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth);
	virtual void reset() override;

	unsigned long getTotalBitSize() const;
	unsigned long getTotalByteSize() const;
	unsigned long getTotalKByteSize() const;
	unsigned long getTotalMByteSize() const;

	unsigned long getTotalLeafBitSize() const;
	unsigned long getTotalLeafByteSize() const;
	double getTotalLeafMByteSize() const;

	unsigned long getTotalLhcBitSize() const;
	unsigned long getTotalLhcByteSize() const;
	unsigned long getTotalLhcKByteSize() const;
	unsigned long getTotalLhcMByteSize() const;

	unsigned long getTotalAhcBitSize() const;
	unsigned long getTotalAhcByteSize() const;
	unsigned long getTotalAhcKByteSize() const;
	unsigned long getTotalAhcMByteSize() const;

protected:
	std::ostream& output(std::ostream &out) const override;

private:
	unsigned long totalLHCByteSize;
	unsigned long totalAHCByteSize;
	unsigned long totalTreeByteSize;
	unsigned long totalLeafByteSize;

	template <unsigned int PREF_BLOCKS>
	unsigned long superSize(const TNode<DIM, PREF_BLOCKS>* node);
};

#include "../visitors/SizeVisitor.h"
#include "../nodes/Node.h"
#include "../nodes/TNode.h"
#include "../nodes/LHC.h"
#include "../nodes/AHC.h"

using namespace std;

template <unsigned int DIM>
SizeVisitor<DIM>::SizeVisitor() : Visitor<DIM>() {
	reset();
}

template <unsigned int DIM>
SizeVisitor<DIM>::~SizeVisitor() { }

template <unsigned int DIM>
template <unsigned int WIDTH>
void SizeVisitor<DIM>::visitSub(PHTree<DIM, WIDTH>* tree) {
	totalTreeByteSize += sizeof (tree->root_);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS, unsigned int N>
void SizeVisitor<DIM>::visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth) {
	totalLHCByteSize += this->template superSize<PREF_BLOCKS>(node);
	totalLHCByteSize += sizeof (node->addresses_);
	totalLHCByteSize += sizeof (node->m);
	totalLHCByteSize += sizeof (node->references_);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void SizeVisitor<DIM>::visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth) {
	totalAHCByteSize += this->template superSize<PREF_BLOCKS>(node);
	totalAHCByteSize += sizeof(node->nContents);
	totalAHCByteSize += sizeof(node->references_);
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
unsigned long SizeVisitor<DIM>::superSize(const TNode<DIM, PREF_BLOCKS>* node) {
	unsigned long prefixSize = sizeof (node->prefixBits_);
	prefixSize += sizeof (node->prefix_);
	unsigned long suffixSize = 0;
	if (node->getSuffixStorage()) {
		suffixSize = node->getSuffixStorage()->getByteSize();
	}
	totalLeafByteSize += suffixSize;

	return prefixSize;
}

template <unsigned int DIM>
void SizeVisitor<DIM>::reset() {
	totalLHCByteSize = 0;
	totalAHCByteSize = 0;
	totalTreeByteSize = 0;
	totalLeafByteSize = 0;
}

template <unsigned int DIM>
std::ostream& SizeVisitor<DIM>::output(std::ostream &out) const {
	float lhcSizePercent = float(totalLHCByteSize) * 100 / float(getTotalByteSize());
	float ahcSizePercent = float(totalAHCByteSize) * 100 / float(getTotalByteSize());
	return out << "total size: " << getTotalKByteSize()
			<< "KByte | " << getTotalMByteSize()
			<< "MByte (LHC: " << lhcSizePercent << "%, AHC: "
			<< ahcSizePercent << "%), suffixes (leafs): " << getTotalLeafMByteSize() << "MByte" << std::endl;
}

template <unsigned int D>
std::ostream& operator <<(std::ostream &out, const SizeVisitor<D>& v) {
	return v.output(out);
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalLeafByteSize() const {
	return totalLeafByteSize;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalLeafBitSize() const {
	return getTotalLeafByteSize() * 8;
}

template <unsigned int DIM>
double SizeVisitor<DIM>::getTotalLeafMByteSize() const {
	return double(getTotalLeafByteSize()) / 1000000.0;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalBitSize() const {
	return getTotalByteSize() * 8;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalByteSize() const {
	return totalLHCByteSize + totalAHCByteSize + totalLeafByteSize;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalKByteSize() const {
	return getTotalByteSize() / 1000;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalMByteSize() const {
	return getTotalByteSize() / 1000000;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalLhcBitSize() const {
	return getTotalLhcByteSize() * 8;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalLhcByteSize() const {
	return totalLHCByteSize;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalLhcKByteSize() const {
	return getTotalLhcByteSize() / 1000;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalLhcMByteSize() const {
	return getTotalLhcByteSize() / 1000000;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalAhcBitSize() const {
	return getTotalAhcByteSize() * 8;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalAhcByteSize() const {
	return totalAHCByteSize;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalAhcKByteSize() const {
	return getTotalAhcByteSize() / 1000;
}

template <unsigned int DIM>
unsigned long SizeVisitor<DIM>::getTotalAhcMByteSize() const {
	return getTotalAhcByteSize() / 1000000;
}

#endif /* SRC_VISITORS_SIZEVISITOR_H_ */
