////////////////////////////////////////////////////////////////////////
// \file    StandardRecord.h
// \brief   StandardRecord defines top-level objects for
//          Common Analysis File trees.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef STANDARDRECORD_H
#define STANDARDRECORD_H

#include "SRCRTHit.h"
#include "SRHeader.h"
#include "SRSlice.h"
#include "SRSliceRecoBranch.h"
#include "SRTruthBranch.h"

/// Common Analysis Files
namespace caf
{

  /// \brief   The StandardRecord is the primary top-level object in the
  ///          Common Analysis File trees.

  class StandardRecord
  {

  public:
    StandardRecord();
    ~StandardRecord();

    SRHeader          hdr;    ///< Header branch: run, subrun, etc.
    //    SRSpill          spill;  ///< Beam spill branch: pot, beam current, etc.
    SRSlice           slc;    ///< Slice branch: nhit, extents, time, etc.
    SRSliceRecoBranch reco;   ///< Slice reco branch: tracks, showers, etc.
    SRTruthBranch     mc;     ///< Slice truth branch for matching interaction

    std::vector<SRTrueParticle> true_particles; ///< List of true particles in the spill
    std::vector<SRCRTHit> crt_hits; ///!< List of CRT hits in the event

  };

} // end namespace

#endif // STANDARDRECORD_H
//////////////////////////////////////////////////////////////////////////////
