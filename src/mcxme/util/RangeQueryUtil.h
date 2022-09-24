/*
 * RangeQueryUtil.h
 *
 *  Created on: May 24, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_RANGEQUERYUTIL_H_
#define SRC_UTIL_RANGEQUERYUTIL_H_

#include <math.h>

#include "../PHTree.h"
#include "../Entry.h"
#include "../iterators/RangeQueryIterator.h"
#include "../util/RandUtil.h"

template <unsigned int DIM, unsigned int WIDTH>
class RangeQueryUtil {
public:

	static RangeQueryIterator<DIM, WIDTH>* getFullRangeIterator(
			const PHTree<DIM, WIDTH>& tree) {
		return getSkewedRangeIterator(tree, 0.0, 1.0);
	}

	static RangeQueryIterator<DIM, WIDTH>* getSelectiveRangeIteratorRandom(
			const PHTree<DIM, WIDTH>& tree,
			double selectivity) {

		assert (selectivity > 0.0 && selectivity <= 1.0);

		const double perDimSelectivity = pow(selectivity, 1.0 / double(DIM));
		assert (perDimSelectivity > selectivity && perDimSelectivity <= 1.0);

		const unsigned long domainSize = (1uL << WIDTH) - 1uL;
		const unsigned long hyperRectSize = domainSize * perDimSelectivity;
		const unsigned long domainWidth = domainSize - hyperRectSize;

		std::vector<unsigned long> lower = RandUtil::generateRandValues(DIM, 0, domainWidth);
		std::vector<unsigned long> upper(DIM);
		for (unsigned d = 0; d < DIM; ++d) {
			upper[d] = lower[d] + hyperRectSize;
			assert (lower[d] < upper[d] && upper[d] < domainSize);
		}

		return tree.rangeQuery(lower, upper);
	}

	static RangeQueryIterator<DIM, WIDTH>* getSkewedRangeIterator(
				const PHTree<DIM, WIDTH>& tree, double fromLowerOffsetPercent, double toUpperOffsetPercent) {
			vector<unsigned long> lower(DIM);
			vector<unsigned long> upper(DIM);

			assert (fromLowerOffsetPercent >= 0.0 && fromLowerOffsetPercent <= 1.0);
			assert (toUpperOffsetPercent >= 0.0 && toUpperOffsetPercent <= 1.0);
			unsigned long maxValue = (WIDTH == 64)? -1uL : (1uL << WIDTH) - 1uL;
			unsigned long lowerValue = floor(double(maxValue) * fromLowerOffsetPercent);
			unsigned long upperValue = ceil(double(maxValue) * toUpperOffsetPercent);
			assert (lowerValue <= upperValue);
			for (unsigned d = 0; d < DIM; ++d) {
				lower.at(d) = lowerValue;
				upper.at(d) = upperValue;
			}

			return tree.rangeQuery(lower, upper);
		}

	static unsigned int countEntriesInFullRange(const PHTree<DIM, WIDTH>& tree) {
		RangeQueryIterator<DIM, WIDTH>* it = getFullRangeIterator(tree);
		unsigned int count = countEntriesInRange(it, tree);
		delete it;
		return count;
	}

	static bool fullRangeContainsId(const PHTree<DIM, WIDTH>& tree, int id) {
		RangeQueryIterator<DIM, WIDTH>* it = getFullRangeIterator(tree);
		bool contained = rangeContainsId(it, tree, id);
		delete it;
		return contained;
	}

	static bool rangeContainsId(const Entry<DIM, WIDTH>& lowerLeft, const Entry<DIM, WIDTH>& upperRight,
			const PHTree<DIM, WIDTH>& tree, int id) {
		RangeQueryIterator<DIM, WIDTH>* it = tree.rangeQuery(lowerLeft, upperRight);
		bool contained = rangeContainsEntry(it, tree, id);
		delete it;
		return contained;
	}

	static bool rangeContainsId(RangeQueryIterator<DIM, WIDTH>* it,
				const PHTree<DIM, WIDTH>& tree, int id) {

			bool foundEqualEntry = false;
			while (it->hasNext() && !foundEqualEntry) {
				Entry<DIM, WIDTH> entryInRange = it->next();
				foundEqualEntry = id == entryInRange.id_;
			}

			return foundEqualEntry;
		}

	static unsigned int countEntriesInRange(const Entry<DIM, WIDTH>& lowerLeft, const Entry<DIM, WIDTH>& upperRight,
				const PHTree<DIM, WIDTH>& tree) {
		RangeQueryIterator<DIM, WIDTH>* it = tree.rangeQuery(lowerLeft, upperRight);
		unsigned int count = countEntriesInRange(it, tree);
		delete it;
		return count;
	}

	static unsigned int countEntriesInRange(RangeQueryIterator<DIM, WIDTH>* it,
			const PHTree<DIM, WIDTH>& tree) {
		unsigned int nEntriesInRange = 0;
		while (it->hasNext()) {
			Entry<DIM, WIDTH> entryInRange = it->next();
			++nEntriesInRange;
		}

		return nEntriesInRange;
	}
};

#endif /* SRC_UTIL_RANGEQUERYUTIL_H_ */
