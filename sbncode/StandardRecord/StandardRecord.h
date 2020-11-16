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
#include "SRFakeReco.h"

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
    SRSliceRecoBranch reco;   ///< Slice reco branch: tracks, showers, etc.
    SRTruthBranch     mc;     ///< Truth branch for all interactions

    int                        nslc;    ///< Number of slices in list
    std::vector<SRSlice>        slc;    ///< Slice branch.
    int                        nfake_reco; ///< Number of Fake-Reco's in list
    std::vector<SRFakeReco>     fake_reco; ///< List of fake-reco slices
    int                        ntrue_particles; ///< Number of true particles in list
    std::vector<SRTrueParticle> true_particles; ///< True particles in spill
    int                        ncrt_hits; ///!< Number of CRT hits in event
    std::vector<SRCRTHit>       crt_hits; ///!< CRT hits in event


    bool pass_flashtrig;     ///< Whether this Record passed the Flash Trigger requirement

  };

} // end namespace

#endif // STANDARDRECORD_H
//////////////////////////////////////////////////////////////////////////////
