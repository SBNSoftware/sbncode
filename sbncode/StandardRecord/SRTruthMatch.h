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

      float      visEinslc;     ///< True deposited energy in slice [GeV]
      float      visEcosmic;    ///< True slice deposited energy from cosmics
      float      eff;           ///< Slice efficiency for this interaction
      float      pur;           ///< Slicer purity for this interaction
      int        index;         ///< Index of the matched true neutrino interaction (-1 if not matched to neutrino)
      bool       is_numucc_primary; ///< Whether this is the "primary" reco neutrino slice as defined by the numu CC analysis

      void setDefault();

    };
} // end namespace

#endif // SRTRUTHMATCH_H
//////////////////////////////////////////////////////////////////////////////
