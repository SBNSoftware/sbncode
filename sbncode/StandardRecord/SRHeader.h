////////////////////////////////////////////////////////////////////////
// \file    SRHeader.h
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRHEADER_H
#define SRHEADER_H

#include "sbncode/StandardRecord/SREnums.h"

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
      unsigned int   fno;       ///< Index of file processed by CAFMaker
      unsigned int   ngenevt;   ///< Number of events generated in input file associated with this record (before any filters)
      float          pot;       ///< POT of the subrun associated with this record
      MCType_t       mctype;    ///< Type of Monte Carlo used to generate this record
      Det_t          det;       ///< Detector, SBND=1 ICARUS=2
      bool           first_in_subrun; ///< Whether this event is the first in the subrun
      bool           first_in_file;   ///< Whether this event is there first in the input file
      // bool           blind;     ///< if true, record has been corrupted for blindness

    };

} // end namespace

#endif // SRHEADER_H
//////////////////////////////////////////////////////////////////////////////
