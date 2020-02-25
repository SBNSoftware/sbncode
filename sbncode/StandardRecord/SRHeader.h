////////////////////////////////////////////////////////////////////////
// \file    SRHeader.h
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRHEADER_H
#define SRHEADER_H

#include "SREnums.h"

namespace caf
{
  /// Header representing overview information for the current event/slice
  class SRHeader
    {
    public:
      SRHeader();
      ~SRHeader();

      unsigned int   run;       ///< run number
      unsigned int   subrun;    ///< subrun number
      // int            cycle;     ///< MC simulation cycle number
      // int            batch;     ///< MC simulation batch number
      unsigned int   evt;       ///< ART event number, indexes trigger windows.
      unsigned short subevt;    ///< slice number within spill
      bool           ismc;      ///< data or MC?  True if MC
      Det_t          det;       ///< Detector, SBND=1 ICARUS=2
      // bool           blind;     ///< if true, record has been corrupted for blindness

    };

} // end namespace

#endif // SRHEADER_H
//////////////////////////////////////////////////////////////////////////////
