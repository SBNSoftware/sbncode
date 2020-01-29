////////////////////////////////////////////////////////////////////////
// \file    StandardRecord.h
// \brief   StandardRecord defines top-level objects for
//          Common Analysis File trees.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef STANDARDRECORD_H
#define STANDARDRECORD_H

#include "SRHeader.h"
#include "SRSlice.h"
#include "SRSliceRecoBranch.h"
#include "SRParticle.h"

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

    SRHeader         hdr;    ///< Header branch: run, subrun, etc.
    //    SRSpill          spill;  ///< Beam spill branch: pot, beam current, etc.
    SRSlice          slc;    ///< Slice branch: nhit, extents, time, etc.
    SRSliceRecoBranch reco;    ///< Slice reco branch: tracks, showers, etc.

    // TODO (from Gray): How do we represent true particles here??
    // The true particles are per-event, I don't see how they fit into
    // a per-slice structure
    std::vector<SRParticle> particle;

  };

} // end namespace

#endif // STANDARDRECORD_H
//////////////////////////////////////////////////////////////////////////////
