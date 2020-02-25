////////////////////////////////////////////////////////////////////////
// \file    SRSliceRecoBranch.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICERECOBRANCH_H
#define SRSLICERECOBRANCH_H

#include "SRTrack.h"
#include "SRShower.h"
// #include "SRVertex.h"

#include <vector>

namespace caf
{
  /// Vectors of reconstructed vertices found by various algorithms
  class SRSliceRecoBranch
  {
  public:
    SRSliceRecoBranch();
    ~SRSliceRecoBranch();
            
    /* std::vector<SRVertex> vtx; ///< Vector of vertices */
    /* size_t               nvtx;   ///< Number of vertices */

    std::vector<SRTrack>  trk;      ///< Vector of pandora tracks
    size_t               ntrk;     ///< Number of panora tracks

    std::vector<SRShower> shw;      ///< Vector of pandora showers
    size_t               nshw;     ///< Number of panora showers

    void fillSizes();
      
  };
  
} // end namespace

#endif // SRSLICERECOBRANCH_H
////////////////////////////////////////////////////////////////////////////
