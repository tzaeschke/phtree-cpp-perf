#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <sys/time.h>
#include <unordered_map>

#include "ctpl_stl.h"
#include "BBTree.h"

/**
 * Determines the current time; may be used for measuring the execution time
 * of a code segment.
 *
 * Example:
 * double start = gettime();
 * // code segment
 * double execution_time = gettime() - start;
 */
static double gettime(void) {
  struct timeval now_tv;
  gettimeofday (&now_tv,NULL);

  return ((double)now_tv.tv_sec) + ((double)now_tv.tv_usec) / 1000000.0;
}

// determines the statistical average (mean) of a given sequence of doubles.
static double getaverage(double* runtimes, size_t n) {
  double sum = 0.0;

  for (size_t i = 0; i < n; ++i) {
    sum += runtimes[i];
  }

  return (sum / n);
}

// determines the standard deviation of a given sequence of doubles.
static double getstddev(double* runtimes, size_t n) {
  double avg = 0.0;
  double std_dev = 0.0;

  for (size_t i = 0; i < n; ++i) {
    avg += runtimes[i];
  }
  avg = (avg / n);

  for (size_t i = 0; i < n; ++i) {
    std_dev += (runtimes[i] - avg) * (runtimes[i] - avg);
  }
  std_dev = std_dev / n;
  std_dev = sqrt(std_dev);

  return std_dev;
}

int main(int argc, char* argv[]) {
  // random seed
  srand(0);
  // default number of data objects
  size_t n = 1000000;
  // default dimensionality of the feature space
  size_t m = 10;
  // default number of to-be-executed range queries
  size_t rq = 100;

  if (argc >= 2)
    n = atoi(argv[1]);
  if (argc >= 3)
    m = atoi(argv[2]);
  std::cout << "n = " << n << ", m = " << m << std::endl;

  // default DOP
  size_t threads = 1;
  // default value range for synthetic data
  int o = 1;
  // default selectivity for range queries
  float selectivity = 0.5;
  if (atoi(argv[3]) == 2 && argc == 5) {
    selectivity = atof(argv[4]);
    std::cout << "Query Selectivity per Dimension: " << selectivity << std::endl;
  }

  std::vector< std::vector<float> > data_points(n, std::vector<float>(m));
  if (argc == 6) {
          threads = atoi(argv[5]);
  } else {
          threads = std::thread::hardware_concurrency();
  }
  std::cout << "THREADS: " << threads << std::endl;

  // initialize
  BBTree bbtree(m, threads);

  // load or generate data objects
  if (atoi(argv[3]) == 3) { // dataset GENOMIC
  } else if(atoi(argv[3]) == 4) { // dataset POWER
  } else if(atoi(argv[3]) == 1) { // dataset CLUSTERED
  } else { // dataset UNIFORM
    for (size_t i = 0; i < n; ++i) {
      std::vector<float> data_point(m);
      for (size_t j = 0; j < m; ++j) {
        data_point[j] = (float) ((rand() % (o * 1000000)) / 1000000.0);
      }
      data_points[i] = data_point;
    }
  }

  std::vector<std::vector<float> > lb_queries(rq,
    std::vector<float>(m, std::numeric_limits<float>::min()));
  std::vector<std::vector<float> > ub_queries(rq,
    std::vector<float>(m, std::numeric_limits<float>::max()));
  // Generate or load rq queries
  if (atoi(argv[3]) == 3) { // dataset GENOMIC
  } else { // generate synthetic range queries
    for (size_t i = 0; i < rq; ++i) {
      int first  = rand() % n;
      int second = rand() % n;
      for (size_t j = 0; j < m; ++j) {
        lb_queries[i][j] = std::min(data_points[first][j],
                                    data_points[second][j]);
        ub_queries[i][j] = std::max(data_points[first][j],
                                    data_points[second][j]);
        if (argc == 5 && atoi(argv[3]) == 2) { // fixed selectivity
          lb_queries[i][j] = (float) ((rand() % (o * 10000)) /
              1000000.0);
          ub_queries[i][j] = lb_queries[i][j] + o * selectivity;
        }
      }
    }
  }

  // shuffle data objects to achieve a random insertion order
  std::random_shuffle(data_points.begin(), data_points.end());

  double *runtimes;
  double start;

  std::cout << "BB-Tree [inserts]" << std::endl;
  runtimes = new double[n];
  for (size_t i = 0; i < n; ++i) {
    start = gettime();
    bbtree.InsertObject(data_points[i], (uint32_t) (i+1));
    runtimes[i] = (gettime() - start) * 1000000;
  }
  std::cout << "Mean: " << getaverage(runtimes, n) << " Standard Deviation: " <<
               getstddev(runtimes, n) << std::endl;
  delete runtimes;

  std::cout << "BB-Tree [point queries]" << std::endl;
  runtimes = new double[n];
  for (size_t i = 0; i < n; ++i) {
    start = gettime();
    assert((i+1) == bbtree.SearchObject(data_points[i]));
    runtimes[i] = (gettime() - start) * 1000000;
  }
  std::cout << "Mean: " << getaverage(runtimes, n) << " Standard Deviation: " <<
               getstddev(runtimes, n) << std::endl;
  delete runtimes;

  std::cout << "BB-Tree [range queries]" << std::endl;
  runtimes = new double[rq];
  for (size_t i = 0; i < rq; ++i) {
    start = gettime();
    std::vector<uint32_t> results =  bbtree.SearchRange(lb_queries[i],
                                                       ub_queries[i]);
    runtimes[i] = (gettime() - start) * 1000;
  }
  std::cout << "Mean: " << getaverage(runtimes, rq) << " Standard Deviation: " <<
               getstddev(runtimes, rq) << std::endl;
  delete runtimes;

  std::cout << "BB-Tree [range queries/multithreaded]" << std::endl;
  runtimes = new double[rq];
  for (size_t i = 0; i < rq; ++i) {
    start = gettime();
    std::vector<uint32_t> results =  bbtree.SearchRangeMT(lb_queries[i],
                                                         ub_queries[i]);
    runtimes[i] = (gettime() - start) * 1000;
  }
  std::cout << "Mean: " << getaverage(runtimes, rq) << " Standard Deviation: " <<
               getstddev(runtimes, rq) << std::endl;
  delete runtimes;

  std::cout << "BB-Tree [deletes]" << std::endl;
  runtimes = new double[n];
  for (size_t i = 0; i < n; ++i) {
    std::vector<float> object_to_delete = data_points[i];
    assert((i+1) == bbtree.SearchObject(object_to_delete));
    size_t old_count = bbtree.getCount();
    start = gettime();
    assert(true == bbtree.DeleteObject(object_to_delete));
    runtimes[i] = (gettime() - start) * 1000000;
    assert(old_count-1 == bbtree.getCount());
  }
  std::cout << "Mean: " << getaverage(runtimes, n) << " Standard Deviation: " <<
               getstddev(runtimes, n) << std::endl << std::endl;
  delete runtimes;

  return 0;
}
