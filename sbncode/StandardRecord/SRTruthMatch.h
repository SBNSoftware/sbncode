////////////////////////////////////////////////////////////////////////
// \file    SRTruthMatch.h
// \brief   SRTruthMatch object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRTRUTHMATCH_H
#define SRTRUTHMATCH_H

#include "SRVector3D.h"


namespace caf
{
  /// An  SRTruthMatch contains overarching information for a slice.
  class SRTruthMatch
    {
    public:

      SRTruthMatch();
      virtual ~SRTruthMatch();

      float      eff;           ///< Slice efficiency for this interaction
      float      pur;           ///< Slicer purity for this interaction
      int        index;         ///< Index of the matched true neutrino interaction (-1 if not matched to neutrino)

      void setDefault();

    };
} // end namespace

#endif // SRTRUTHMATCH_H
//////////////////////////////////////////////////////////////////////////////
