/*
 * Visitor.h
 *
 *  Created on: Mar 10, 2016
 *      Author: max
 */

#ifndef VISITORS_VISITOR_H_
#define VISITORS_VISITOR_H_

#include <iostream>

template <unsigned int DIM, unsigned int PREF_BLOCKS, unsigned int N>
class LHC;

template <unsigned int DIM, unsigned int PREF_BLOCKS>
class AHC;

template <unsigned int DIM, unsigned int WIDTH>
class PHTree;

template <unsigned int DIM>
class Visitor {
public:
	Visitor();
	virtual ~Visitor();

	template <unsigned int WIDTH>
	void visit(PHTree<DIM, WIDTH>* tree);
	template <unsigned int PREF_BLOCKS>
	void visit(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth, unsigned int index);
	template <unsigned int PREF_BLOCKS, unsigned int N>
	void visit(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth, unsigned int index);
	virtual void reset() =0;
	std::ostream& operator <<(std::ostream &out);

protected:
	virtual std::ostream& output(std::ostream &out) const =0;
};

#include <stdexcept>
#include "../visitors/SizeVisitor.h"
#include "../visitors/PrefixSharingVisitor.h"
#include "../visitors/CountNodeTypesVisitor.h"
#include "../visitors/AssertionVisitor.h"
#include "../visitors/SuffixVisitor.h"

template <unsigned int DIM>
Visitor<DIM>::Visitor() {}

template <unsigned int DIM>
Visitor<DIM>::~Visitor() {}

template <unsigned int DIM>
template <unsigned int WIDTH>
void Visitor<DIM>::visit(PHTree<DIM, WIDTH>* tree) {
	if (SizeVisitor<DIM>* sub = dynamic_cast<SizeVisitor<DIM>*>(this)) {
		sub->visitSub(tree);
	} else if (PrefixSharingVisitor<DIM>* sub = dynamic_cast<PrefixSharingVisitor<DIM>*>(this)) {
		sub->visitSub(tree);
	} else if (CountNodeTypesVisitor<DIM>* sub = dynamic_cast<CountNodeTypesVisitor<DIM>*>(this)) {
		sub->visitSub(tree);
	} else if (AssertionVisitor<DIM>* sub = dynamic_cast<AssertionVisitor<DIM>*>(this)) {
		sub->visitSub(tree);
	} else if (SuffixVisitor<DIM>* sub = dynamic_cast<SuffixVisitor<DIM>*>(this)) {
		sub->visitSub(tree);
	} else {
		throw std::runtime_error("unknown visitor");
	}
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS>
void Visitor<DIM>::visit(AHC<DIM, PREF_BLOCKS>* node, unsigned int depth, unsigned int index) {
	if (SizeVisitor<DIM>* sub = dynamic_cast<SizeVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (PrefixSharingVisitor<DIM>* sub = dynamic_cast<PrefixSharingVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (CountNodeTypesVisitor<DIM>* sub = dynamic_cast<CountNodeTypesVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (AssertionVisitor<DIM>* sub = dynamic_cast<AssertionVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (SuffixVisitor<DIM>* sub = dynamic_cast<SuffixVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth, index);
	} else {
		throw std::runtime_error("unknown visitor");
	}
}

template <unsigned int DIM>
template <unsigned int PREF_BLOCKS, unsigned int N>
void Visitor<DIM>::visit(LHC<DIM, PREF_BLOCKS, N>* node, unsigned int depth, unsigned int index) {
	if (SizeVisitor<DIM>* sub = dynamic_cast<SizeVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (PrefixSharingVisitor<DIM>* sub = dynamic_cast<PrefixSharingVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (CountNodeTypesVisitor<DIM>* sub = dynamic_cast<CountNodeTypesVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (AssertionVisitor<DIM>* sub = dynamic_cast<AssertionVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth);
	} else if (SuffixVisitor<DIM>* sub = dynamic_cast<SuffixVisitor<DIM>*>(this)) {
		sub->visitSub(node, depth, index);
	} else {
		throw std::runtime_error("unknown visitor");
	}
}


template <unsigned int DIM>
std::ostream& Visitor<DIM>::operator <<(std::ostream &out) {
	return output(out);
}


#endif /* VISITORS_VISITOR_H_ */
