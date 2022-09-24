/*
 * FileInputUtil.h
 *
 *  Created on: Mar 15, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_FILEINPUTUTIL_H_
#define SRC_UTIL_FILEINPUTUTIL_H_

#include <vector>

class FileInputUtil {
public:
	// parses the file at the given location in the format 'ulong, ulong, ulong, ...\n...'
	template <unsigned int DIM>
	static std::vector<vector<unsigned long>>* readEntries(std::string fileLocation);
	template <unsigned int DIM>
	static std::vector<vector<double>>* readEntriesToFloat(std::string fileLocation);
	// parses the file at the given location in the format 'float, float, float, ...\n...'
	template <unsigned int DIM>
	static std::vector<vector<unsigned long>>* readFloatEntries(string fileLocation, size_t decimals);
};

#include <iostream>
#include <fstream>
#include <regex>
#include <vector>
#include <limits.h>
#include <string>
#include <assert.h>
#include <stdexcept>
#include <math.h>
#include "../Entry.h"
#include "../util/FileInputUtil.h"

using namespace std;

inline vector<double> getNextLineTokensFloat(ifstream& stream) {
	string line;
	getline(stream, line);
	stringstream lineStream(line);
	string cell;
	vector<double> tokens;

	while (getline(lineStream, cell, ' ')) {
		const double parsedToken = stod(cell);
		tokens.push_back(parsedToken);
	}

	return tokens;
}

inline vector<unsigned long> getNextLineTokens(ifstream& stream) {
	string line;
	getline(stream, line);
	stringstream lineStream(line);
	string cell;
	vector<unsigned long> tokens;

	while (getline(lineStream, cell, ',')) {
		const unsigned long parsedToken = stoi(cell);
		tokens.push_back(parsedToken);
	}

	return tokens;
}

inline vector<unsigned long> getNextLineTokens(ifstream& stream, unsigned long decimals) {
	string line;
	getline(stream, line);
	stringstream lineStream(line);
	string cell;
	vector<unsigned long> tokens;
	const long double decimalShift = pow10(decimals);

	while (getline(lineStream, cell, ' ')) {
		// TODO added a shift
/*		const long double parsedToken = stold(cell) + 2000.0L;
		const long double longToken = parsedToken * decimalShift;
		if (parsedToken < 0.0)
			throw runtime_error("Can only handle positive values: the offset is too small");
		if (parsedToken > nextafter(ULONG_MAX, 0))
			throw runtime_error("Overflow: The floating point value cannot be represented as a 64-bit integer.");
		const unsigned long convertedToken = (unsigned long) longToken;*/
		double value = stod(cell);
		if (value == -0.0) {
			value = 0.0;
		}
		unsigned long convertedToken;
		memcpy(&convertedToken, &value, sizeof(value));
		if (value < 0.0) {
			convertedToken = (~convertedToken) | (1L << 63);
		}

		tokens.push_back(convertedToken);
	}

	return tokens;
}

template <unsigned int DIM>
vector<vector<unsigned long>>* FileInputUtil::readEntries(string fileLocation) {

	ifstream myfile (fileLocation);
	vector<vector<unsigned long>>* result = new vector<vector<unsigned long>>();
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			vector<unsigned long> values = getNextLineTokens(myfile);
			if (!values.empty()) {
				assert (values.size() == DIM);
				result->push_back(values);
			}
		}
	} else {
		delete result;
		throw runtime_error("cannot open the file " + fileLocation);
	}

	return result;
}

template <unsigned int DIM>
vector<vector<double>>* FileInputUtil::readEntriesToFloat(string fileLocation) {

	ifstream myfile (fileLocation);
	vector<vector<double>>* result = new vector<vector<double>>();
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			vector<double> values = getNextLineTokensFloat(myfile);
			if (!values.empty()) {
				assert (values.size() == DIM);
				result->push_back(values);
			}
		}
	} else {
		delete result;
		throw runtime_error("cannot open the file " + fileLocation);
	}

	return result;
}

template <unsigned int DIM>
vector<vector<unsigned long>>* FileInputUtil::readFloatEntries(string fileLocation, size_t decimals) {

	ifstream myfile (fileLocation);
	vector<vector<unsigned long>>* result = new vector<vector<unsigned long>>();
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			vector<unsigned long> values = getNextLineTokens(myfile, decimals);
			if (!values.empty()) {
				assert (values.size() == DIM);
				result->push_back(values);
			}
		}
	} else {
		delete result;
		throw runtime_error("cannot open the file " + fileLocation);
	}

	return result;
}

#endif /* SRC_UTIL_FILEINPUTUTIL_H_ */
