/*
 * RTreeBulkWrapper.h
 *
 *  Created on: Aug 22, 2016
 *      Author: max
 */

#ifndef SRC_UTIL_COMPARE_RTREEBULKWRAPPER_H_
#define SRC_UTIL_COMPARE_RTREEBULKWRAPPER_H_

#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = bg::index;
// 3D - spatial objects
typedef bg::model::point<double, 3, bg::cs::cartesian> point3D;
typedef bg::model::box<point3D> box3D;
typedef std::pair<box3D, size_t> idBox3D;
// 2D - spatial objects
typedef bg::model::point<double, 2, bg::cs::cartesian> point2D;
typedef bg::model::box<point2D> box2D;
typedef std::pair<box2D, size_t> idBox2D;

class RTreeBulkWrapper {
public:
	RTreeBulkWrapper();
	~RTreeBulkWrapper();

	void bulkLoad3DSpatialObject(std::vector<std::vector<double>>& entries);
	void bulkLoad2DSpatialObject(std::vector<std::vector<double>>& entries);
private:
	static const size_t NODE_CAPACITY = 16;

	bgi::rtree<idBox3D, bgi::linear<NODE_CAPACITY> >* rtree3D;
	bgi::rtree<idBox2D, bgi::linear<NODE_CAPACITY> >* rtree2D;
};

#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <vector>
#include <iostream>

using namespace std;



RTreeBulkWrapper::RTreeBulkWrapper() : rtree3D(NULL), rtree2D(NULL) {}

RTreeBulkWrapper::~RTreeBulkWrapper() {
	if (rtree3D) {
		delete rtree3D;
	}

	if (rtree2D) {
		delete rtree2D;
	}
}

void RTreeBulkWrapper::bulkLoad2DSpatialObject(std::vector<std::vector<double>>& entries) {
	assert (entries.size() > 0);
	assert (entries[0].size() == 4);

	// transform to boost data types:
    std::vector<idBox2D> values;
    values.reserve(entries.size());
    for (unsigned e = 0; e < entries.size(); ++e) {
    	point2D lower(entries[e][0], entries[e][1]);
    	point2D upper(entries[e][2], entries[e][3]);
    	box2D b(lower, upper);
    	values.push_back(make_pair(b, e));
    }

    rtree2D = new bgi::rtree<idBox2D, bgi::linear<NODE_CAPACITY> > (values.begin(), values.end());
}

void RTreeBulkWrapper::bulkLoad3DSpatialObject(std::vector<std::vector<double>>& entries) {
	assert (entries.size() > 0);
	assert (entries[0].size() == 6);

	// transform to boost data types:
    std::vector<idBox3D> values;
    values.reserve(entries.size());
    for (unsigned e = 0; e < entries.size(); ++e) {
    	point3D lower(entries[e][0], entries[e][1], entries[e][2]);
    	point3D upper(entries[e][3], entries[e][4], entries[e][5]);
    	box3D b(lower, upper);
    	values.push_back(make_pair(b, e));
    }

    rtree3D = new bgi::rtree<idBox3D, bgi::linear<NODE_CAPACITY> > (values.begin(), values.end());
}

#endif /* SRC_UTIL_COMPARE_RTREEBULKWRAPPER_H_ */
