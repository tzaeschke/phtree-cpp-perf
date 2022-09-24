#include <iostream>
#include <string>
#include <vector>
#include <assert.h>

#ifndef BOOST_THREAD_VERSION
	// requires upwards conversions of shared mutex
#define BOOST_THREAD_VERSION 3
#endif

#include "Entry.h"
#include "PHTree.h"
#include "util/PlotUtil.h"
#include "util/rdtsc.h"
#include "visitors/CountNodeTypesVisitor.h"
#include "iterators/RangeQueryIterator.h"

using namespace std;

int mainSimpleExample() {
	const unsigned int bitLength = 8;
	Entry<2, bitLength> e1({ 74, 21 }, 1);
	Entry<2, bitLength> e2({ 75, 28 }, 2);
	Entry<2, bitLength> e3({ 124, 7 }, 3);
	Entry<2, bitLength> e4({ 65, 19 }, 4);
	Entry<2, bitLength> e5({ 75, 21 }, 5);

	cout << "example (points): 2D, bit length: " << bitLength << endl;
	CountNodeTypesVisitor<2>* visitor = new CountNodeTypesVisitor<2>();
	uint64_t sta = RDTSC();
	PHTree<2, bitLength>* phtree = new PHTree<2, bitLength>();
	phtree->insert(e1);
	cout << "CPU cycles per insert: " << RDTSC() - sta << endl;
	cout << *phtree;
	phtree->accept(visitor);
	cout << *visitor << endl;

	sta = RDTSC();
	phtree->insert(e2);
	cout << "CPU cycles per insert: " << RDTSC() - sta << endl;
	assert (phtree->lookup(e1).second == 1);
	assert (phtree->lookup(e2).second == 2);
	assert (!phtree->lookup(e3).first);
	cout << *phtree;
	visitor->reset();
	phtree->accept(visitor);
	cout << *visitor << endl;

	sta = RDTSC();
	phtree->insert(e3);
	cout << "CPU cycles per insert: " << RDTSC() - sta << endl;
	assert (phtree->lookup(e1).second == 1);
	assert (phtree->lookup(e2).second == 2);
	assert (phtree->lookup(e3).second == 3);
	cout << *phtree;
	visitor->reset();
	phtree->accept(visitor);
	cout << *visitor << endl;

	sta = RDTSC();
	phtree->insert(e4);
	cout << "CPU cycles per insert: " << RDTSC() - sta << endl;
	assert (phtree->lookup(e1).second == 1);
	assert (phtree->lookup(e2).second == 2);
	assert (phtree->lookup(e3).second == 3);
	assert (phtree->lookup(e4).second == 4);
	cout << *phtree;
	visitor->reset();
	phtree->accept(visitor);
	cout << *visitor << endl;

	sta = RDTSC();
	phtree->insert(e5);
	cout << "CPU cycles per insert: " << RDTSC() - sta << endl;
	assert (phtree->lookup(e5).second == 5);
	cout << *phtree;
	visitor->reset();
	phtree->accept(visitor);
	cout << *visitor << endl;

	cout << "The following entries are in the range (0,0) - (100,100):" << endl;
	RangeQueryIterator<2, bitLength>* it = phtree->rangeQuery({0,0}, {100,100});
	while (it->hasNext()) {
		Entry<2, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following entries are in the range (65,10) - (150,20):" << endl;
	it = phtree->rangeQuery({65,10}, {150,20});
	while (it->hasNext()) {
		Entry<2, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following entries are in the range (70,5) - (200,25):" << endl;
	it = phtree->rangeQuery({70,5}, {200,25});
	while (it->hasNext()) {
		Entry<2, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following entries are in the range (74,20) - (75,21):" << endl;
	it = phtree->rangeQuery({74,20}, {75,21});
	while (it->hasNext()) {
		Entry<2, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	delete visitor;
	delete phtree;

	return 0;
}

void mainBulkExample() {
	const unsigned int bitLength = 6;
	PHTree<1, bitLength>* phtree = new PHTree<1, bitLength>();
	vector<vector<unsigned long>>* values = new vector<vector<unsigned long>>();
	vector<int>* ids = new vector<int>();
	values->push_back({15});	ids->push_back(1);
	values->push_back({18});	ids->push_back(2);
	values->push_back({9});		ids->push_back(3);
	values->push_back({21});	ids->push_back(4);
	values->push_back({14});	ids->push_back(5);
	values->push_back({8}); 	ids->push_back(6);
	values->push_back({19});	ids->push_back(7);
	phtree->bulkInsert(*values, *ids);

	cout << *phtree << endl;

	delete phtree;
	delete values;
	delete ids;
}

int mainFull1DExample() {
	const unsigned int bitLength = 4;
	PHTree<1, bitLength>* phtree = new PHTree<1, bitLength>();

	const unsigned long upperBoundary = (1uL << bitLength);
	for (unsigned long i = 0; i < upperBoundary; ++i) {
		phtree->insert({i}, i);
		assert (phtree->lookup({i}).first);
	}

	cout << *phtree;

	for (unsigned long width = 1; width < upperBoundary; ++width) {
		for (unsigned long lower = 0; lower < upperBoundary - width; ++lower) {
			RangeQueryIterator<1, bitLength>* it = phtree->rangeQuery({lower}, {lower + width});
			unsigned int nIntersects = 0;
			while (it->hasNext()) {
				it->next();
				nIntersects++;
			}

			assert (nIntersects == (1 + width));
			delete it;
		}
	}

	delete phtree;
	return 0;
}

int mainSharing1DExample() {
	const unsigned int bitLength = 6;
	PHTree<1, bitLength>* phtree = new PHTree<1, bitLength>();

	vector<unsigned int> values = {2, 6, 18, 22, 34, 38, 50, 54};
	for (unsigned i = 0; i < values.size(); ++i) {
		const unsigned int value = values[i];
		phtree->insert({value}, value);
		phtree->insert({value + 1}, value + 1);
	}

	cout << (*phtree) << endl;

	RangeQueryIterator<1, bitLength>* it = phtree->rangeQuery({0}, {63});
	unsigned int points = 0;
	while (it->hasNext()) {
		it->next();
		points++;
	}
	assert (points == 16);
	delete it;

	it = phtree->rangeQuery({18}, {39});
	points = 0;
	while (it->hasNext()) {
		it->next();
		points++;
	}
	assert (points == 8);
	delete it;

	delete phtree;
	return 1;
}

int mainHyperCubeExample() {
	const unsigned int bitLength = 4;
	vector<unsigned long> e1Lower = {5, 5};
	vector<unsigned long> e1Upper = {15, 10};
	vector<unsigned long> e2Lower = {5, 10};
	vector<unsigned long> e2Upper = {10, 15};
	vector<unsigned long> e3Lower = {1, 9};
	vector<unsigned long> e3Upper = {6, 14};
	vector<unsigned long> e4Lower = {0, 0};
	vector<unsigned long> e4Upper = {5, 5};

	cout << "example (rectangles): 2D, bit length: " << bitLength << endl;
	PHTree<4, bitLength>* phtree = new PHTree<4, bitLength>();
	phtree->insertHyperRect(e1Lower, e1Upper, 1);
	assert (phtree->lookupHyperRect(e1Lower, e1Upper).second == 1);
	phtree->insertHyperRect(e2Lower, e2Upper, 2);
	assert (phtree->lookupHyperRect(e2Lower, e2Upper).second == 2);
	phtree->insertHyperRect(e3Lower, e3Upper, 3);
	assert (phtree->lookupHyperRect(e3Lower, e3Upper).second == 3);
	phtree->insertHyperRect(e4Lower, e4Upper, 4);
	assert (phtree->lookupHyperRect(e4Lower, e4Upper).second == 4);
	cout << *phtree;

	cout << "The following rectangles touch the range (1,1) - (15,15):" << endl;
	RangeQueryIterator<4, bitLength>* it = phtree->intersectionQuery({1, 1}, {15, 15});
	while (it->hasNext()) {
		Entry<4, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following rectangles touch the range (5,9) - (6,10):" << endl;
	it = phtree->intersectionQuery({5, 9}, {6, 10});
	while (it->hasNext()) {
		Entry<4, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following rectangles touch the range (11,2) - (12,14):" << endl;
	it = phtree->intersectionQuery({11, 2}, {12, 14});
	while (it->hasNext()) {
		Entry<4, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following rectangles are completely inside the range (1,1) - (15,15):" << endl;
	it = phtree->inclusionQuery({1, 1}, {15, 15});
	while (it->hasNext()) {
		Entry<4, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following rectangles are completely inside the range (5,10) - (10,15):" << endl;
	it = phtree->inclusionQuery({5, 10}, {10, 15});
	while (it->hasNext()) {
		Entry<4, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}
	delete it;

	cout << "The following rectangles are completely inside the range (11,2) - (12,14):" << endl;
	it = phtree->inclusionQuery({11, 2}, {12, 14});
	while (it->hasNext()) {
		Entry<4, bitLength> entryInRange = it->next();
		cout << entryInRange << endl;
	}

	delete it;
	delete phtree;
	return 0;
}

int main(int argc, char* argv[]) {

	string debug = "debug";
	string plot = "plot";
	string rand = "rand";
	string benchmark = "benchmark";
	string axon = "axon";

	#ifndef NDEBUG
		cout << "assertions enabled!" << endl;
	#endif

	#ifdef PRINT
		cout << "printing enabled!" << endl;
	#endif

	if (argc != 2 || debug.compare(argv[1]) == 0) {
		mainFull1DExample();
		cout << endl;
		mainSharing1DExample();
		mainSimpleExample();
		cout << endl;
		mainHyperCubeExample();
		cout << endl;
		mainBulkExample();
		return 0;
	} else if (plot.compare(argv[1]) == 0) {
		PlotUtil::plotAverageInsertTimePerDimension<3,64>("./random-extract.dat", false, false);
		PlotUtil::plotAverageInsertTimePerDimension<4,64>("./CA_streets-out-extract.dat", false, true);
		PlotUtil::plotAverageInsertTimePerDimension<6,64>("./axons-extract.dat", false, true);
		PlotUtil::plotAverageInsertTimePerDimension<9,64>("./ped09-out-extract.dat", false, true);
//		PlotUtil::plotAverageInsertTimePerDimension<16,64>("./rea16-out-extract.dat", false, true);

		PlotUtil::plotParallelInsertPerformance<3,64>("./random-extract.dat", false);
		PlotUtil::plotParallelInsertPerformance<4,64>("./CA_streets-out-extract.dat", true);
		PlotUtil::plotParallelInsertPerformance<6,64>("./axons-extract.dat", true);
		PlotUtil::plotParallelInsertPerformance<9,64>("./ped09-out-extract.dat", true);
		PlotUtil::plotParallelInsertPerformance<16,64>("./rea16-out-extract.dat", true);

//		PlotUtil::plotCompareToRTreeBulk<6,64>("./axons.dat", true);
//		PlotUtil::plotCompareParallelTreeToScanQuery<6,64>("./axons.dat", "./ranges.dat", true);
//		PlotUtil::plotParallelInsertPerformance<6,64>("/media/max/TOSHIBA/MA/data/ph-tree_workload/100K-axon-mbr-644000.txt", true);
//		PlotUtil::plotParallelInsertPerformance<3,32>("./benchmark_Java-extract_1M_3D_32bit.dat", false);
//		PlotUtil::plotInsertPerformanceDifferentOrder<6, 64>("./axons.dat", true);
//		PlotUtil::plotInsertPerformanceDifferentOrder<3, 32>("./benchmark_Java-extract_1M_3D_32bit.dat", false);
//		PlotUtil::plotTimeSeriesOfInserts();
//		PlotUtil::plotAverageInsertTimePerDimensionRandom();
//		PlotUtil::plotAverageInsertTimePerNumberOfEntriesRandom();
//		PlotUtil::plotRangeQueryTimePerPercentFilledRandom();
//		PlotUtil::plotRangeQueryTimePerSelectivityRandom();
//		PlotUtil::plotAverageInsertTimePerNumberOfEntries<6, 64>("./axons.dat", true);
		return 0;
	} else if (rand.compare(argv[1]) == 0) {
		vector<size_t> nEntries;
		nEntries.push_back(500000);
		PlotUtil::plotAverageInsertTimePerNumberOfEntriesRandom(nEntries);
	} else if (benchmark.compare(argv[1]) == 0) {
		cout << "run a benchmark extracted from the Java implementation with 1M 3D 32-bit entries" << endl;
		PlotUtil::plotAverageInsertTimePerDimension<3, 32>("./benchmark_Java-extract_1M_3D_32bit.dat", true);
	} else if (axon.compare(argv[1]) == 0) {
		vector<string> axonFiles;
		axonFiles.push_back("./axons.dat");
//		axonFiles.push_back("./dendrites.dat");
		vector<string> dendriteFiles;
		dendriteFiles.push_back("./dendrites.dat");
//		dendriteFiles.push_back("./axons.dat");
		PlotUtil::plotAxonsAndDendrites<6, 64>(axonFiles, dendriteFiles, true);
	} else {
		cerr << "Missing command line argument!" << endl << "valid: 'debug', 'plot', 'rand', 'benchmark', 'axon";
		return 1;
	}
};
