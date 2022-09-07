/**
 * \file seaDBSCAN.h
 *
 * 
 * \brief Class def header for a class seaDBSCAN
 *
 * @author mark ross-lonergan markrl@nevis.columbia.edu
 * Written 20th May 2019.
 */

#ifndef SBNCODE_SINGLEPHOTONANALYSIS_seaDBSCAN_H
#define SBNCODE_SINGLEPHOTONANALYSIS_seaDBSCAN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <climits>
#include <limits>


namespace seaview{
class seaDBSCAN{

  public:
    double m_eps;
    int m_minpts;

    /// constructor
    seaDBSCAN(double in_eps, int in_minpts): m_eps(in_eps), m_minpts(in_minpts) {}

    /// Default destructor
    ~seaDBSCAN(){}

    // identify each hit to a certain cluster, or identify as noise
    std::vector<int> Scan2D(std::vector<std::vector<double>> &pts);
    // get neighbour points for point i from vector pts
    std::vector<std::vector<double>> GetNeighbours(size_t i, std::vector<std::vector<double>> &pts,bool);
    // combine elements in pts to seed if not found in seed
    int UnionSets(std::vector<std::vector<double>> &seed, std::vector<std::vector<double>> &pts);
    // calculate distance between (w1, t1) and (w2, t2)
    double SimpleDist(double w1, double t1, double w2, double t2);

};

}

#endif // SBNCODE_SINGLEPHOTONANALYSIS_seaDBSCAN_H
