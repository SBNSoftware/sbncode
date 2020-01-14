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
//#include "SRSpill.h"
// #include "SRTrackBranch.h"

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
    //    SRTrackBranch    trk;    ///< Track branch: nhit, len, etc.
  };

} // end namespace

#endif // STANDARDRECORD_H
//////////////////////////////////////////////////////////////////////////////
