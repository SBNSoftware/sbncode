////////////////////////////////////////////////////////////////////////
// \file    SRSliceRecoBranch.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICERECOBRANCH_H
#define SRSLICERECOBRANCH_H

#include "sbncode/StandardRecord/SRTrack.h"
#include "sbncode/StandardRecord/SRShower.h"
// #include "sbncode/StandardRecord/SRVertex.h"

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

    std::vector<SRShower> shw;      ///< Vector of trac showers
    size_t               nshw;     ///< Number of trac showers

    // std::vector<SRShower> shw_em;      ///< Vector of em showers
    // size_t               nshw_em;     ///< Number of em showers

    // std::vector<SRShower> shw_pandora;      ///< Vector of pandora showers
    // size_t               nshw_pandora;     ///< Number of pandora showers

    void fillSizes();
      
  };
  
} // end namespace

#endif // SRSLICERECOBRANCH_H
////////////////////////////////////////////////////////////////////////////
