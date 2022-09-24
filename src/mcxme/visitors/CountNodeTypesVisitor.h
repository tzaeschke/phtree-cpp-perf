/*
 * CountNodeTypesVisitor.h
 *
 *  Created on: Mar 10, 2016
 *      Author: max
 */

#ifndef VISITORS_COUNTNODETYPESVISITOR_H_
#define VISITORS_COUNTNODETYPESVISITOR_H_

#include <iostream>
#include "../visitors/Visitor.h"

template <unsigned int DIM>
class CountNodeTypesVisitor: public Visitor<DIM> {
	template <unsigned int D>
	friend std::ostream& operator <<(std::ostream &out, const CountNodeTypesVisitor<D>& v);
public:
	CountNodeTypesVisitor();
	virtual ~CountNodeTypesVisitor();

	template <unsigned int WIDTH>
	void visitSub(PHTree<DIM, WIDTH>* tree);
	template <unsigned int PREF_BLOCKS, unsigned int N>
	void visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth);
	template <unsigned int PREF_BLOCKS>
	void visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth);
	virtual void reset() override;
	std::ostream& operator <<(std::ostream &out) const;

	unsigned long getNumberOfVisitedAHCNodes() const;
	unsigned long getNumberOfVisitedLHCNodes() const;

protected:
	std::ostream& output(std::ostream &out) const override;

private:
	unsigned long nAHCNodes_;
	unsigned long nLHCNodes_;
	vector<unsigned int> lhcSizeHistogram;
};

using namespace std;

/* implementation of template */
template <unsigned int DIM>
CountNodeTypesVisitor<DIM>::CountNodeTypesVisitor() : Visitor<DIM>() {
	reset();
}

template <unsigned int DIM>
CountNodeTypesVisitor<DIM>::~CountNodeTypesVisitor() {
}

template <unsigned int DIM>
void CountNodeTypesVisitor<DIM>::reset() {
	nAHCNodes_ = 0;
	nLHCNodes_ = 0;
	lhcSizeHistogram.clear();
}

template <unsigned int DIM>
template <unsigned int WIDTH>
void CountNodeTypesVisitor<DIM>::visitSub(PHTree<DIM, WIDTH>* tree) {}

template <unsigned int DIM>
unsigned long CountNodeTypesVisitor<DIM>::getNumberOfVisitedAHCNodes() const {
	return nAHCNodes_;
}

template <unsigned int DIM>
unsigned long CountNodeTypesVisitor<DIM>::getNumberOfVisitedLHCNodes() const {
	return nLHCNodes_;
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS, unsigned int N>
void CountNodeTypesVisitor<DIM>::visitSub(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth) {
	nLHCNodes_++;

	while (N >= lhcSizeHistogram.size()) {
		lhcSizeHistogram.push_back(0);
	}

	++lhcSizeHistogram[N];
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void CountNodeTypesVisitor<DIM>::visitSub(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth) {
	nAHCNodes_++;
}

template <unsigned int D>
std::ostream& operator <<(std::ostream &out, const CountNodeTypesVisitor<D>& v) {
	return v.output(out);
}

template <unsigned int DIM>
std::ostream& CountNodeTypesVisitor<DIM>::output(std::ostream &out) const {
	out << "nodes: " << (getNumberOfVisitedAHCNodes() + getNumberOfVisitedLHCNodes());
	out << " (AHC nodes: " << getNumberOfVisitedAHCNodes();
	out << " | LHC nodes: " << getNumberOfVisitedLHCNodes()  << ")"<< endl;
	out << "LHC size histogram:" << endl;
	for (unsigned i = 0; i < lhcSizeHistogram.size(); ++i) {
		out << "\t" << i << ": " << lhcSizeHistogram[i] << endl;
	}
	return out;
}

#endif /* VISITORS_COUNTNODETYPESVISITOR_H_ */
