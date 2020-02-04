////////////////////////////////////////////////////////////////////////
// \file    SRSlice.h
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICE_H
#define SRSLICE_H

#include "SRTrueInteraction.h"
#include <TVector3.h>


namespace caf
{
  /// An  SRSlice contains overarching information for a slice.
  class SRSlice
    {
    public:

      /// A matching of TPC slice charge to Optical flash light
      struct FlashMatch {
        float score;
        float time;
      };


      SRSlice();
      virtual ~SRSlice();


      int      id;          ///< number of hits
      float    charge;      ///< Calorimetric energy
      TVector3 vertex;      ///< Candidate neutrino vertex in local detector coordinates [cm]

      SRTrueInteraction truth; //!< Truth information on the slice      

      void setDefault();

      // To Do: Make SRFlashMatch and add here
      FlashMatch fmatch; //!< Optical flash-match for this slice of TPC charge

      // To Do: Move to slice reco branch 
      bool is_clear_cosmic; //!< Whether pandora marks the slice as a "clear" cosmic
      float nu_score; //!< Score of how neutrino-like the slice is 

    };
} // end namespace

#endif // SRSLICE_H
//////////////////////////////////////////////////////////////////////////////
