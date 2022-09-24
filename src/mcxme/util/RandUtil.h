/*
 * RandUtil.h
 *
 *  Created on: Jun 1, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_RANDUTIL_H_
#define SRC_UTIL_RANDUTIL_H_

#include <math.h>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include <vector>
#include <assert.h>
#include <iostream>
#include <random>

class RandUtil {
public:

	static void init() {
		srand(time(NULL));
	}

	static std::vector<unsigned long> generateRandValues(const size_t nValues) {
		std::vector<unsigned long> randValues(nValues);
		for (unsigned i = 0; i < nValues; ++i) {
			randValues[i] = generateRandValue();
		}

		return randValues;
	}

	static std::vector<unsigned long> generateRandValues(const size_t nValues,
			unsigned long min, unsigned long max) {
		std::vector<unsigned long> randValues(nValues);
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 eng(rd()); // seed the generator
		std::uniform_int_distribution<long> distr(min, max); // define the range
		for (unsigned i = 0; i < nValues; ++i) {
			randValues[i] = static_cast<unsigned long>(distr(eng));
			assert(min <= randValues[i] && randValues[i] <= max);
		}

		return randValues;
	}

	static unsigned long generateRandValue() {
		if (sizeof(unsigned int) < sizeof(unsigned long))
			return (static_cast<unsigned long>(rand())
					<< (sizeof(unsigned int) * 8)) | rand();

		return rand();
	}

	static unsigned long generateRandValueInRange(
			unsigned long min,
			unsigned long max) {
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 eng(rd()); // seed the generator
		std::uniform_int_distribution<long> distr(min, max); // define the range
		long r = distr(eng);
		return static_cast<unsigned long>(r);
	}
};

#endif /* SRC_UTIL_RANDUTIL_H_ */
