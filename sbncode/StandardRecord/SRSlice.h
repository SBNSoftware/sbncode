////////////////////////////////////////////////////////////////////////
// \file    SRSlice.h
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICE_H
#define SRSLICE_H

namespace caf
{
  /// An  SRSlice contains overarching information for a slice.
  class SRSlice
    {
    public:
      SRSlice();
      virtual ~SRSlice();


      int    id;            ///< number of hits

      float  charge;             ///< Calorimetric energy

      void setDefault();

    };
} // end namespace

#endif // SRSLICE_H
//////////////////////////////////////////////////////////////////////////////
