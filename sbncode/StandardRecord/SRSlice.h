////////////////////////////////////////////////////////////////////////
// \file    SRSlice.h
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICE_H
#define SRSLICE_H
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


      /// Truth information on a slice of TPC charge
      struct Truth {
        interaction_mode_ mode; //!< True mode of the TPC charge
        int interaction_id; //!< ID of the truth-matched neutrino interaction (-1 if no such match)
        float total_deposited_energy; //!< Total true deposited energy associated with this slice [GeV]
        float neutrino_match_energy; //!< Totral true deposited energy in this slice associated with the
                                     // neutrino match (0 if no such match) [GeV]
        float cosmic_energy; //!< Total true deposited energy associated with cosmic activity
        float unmatched_energy; //!< Total true deposited energy un-matched with any true particle activity
      };

      SRSlice();
      virtual ~SRSlice();


      int    id;            ///< number of hits

      float  charge;             ///< Calorimetric energy

      TVector3 vertex; //!< Candidate neutrino vertex in local detector coordinates [cm]
      FlashMatch fmatch; //!< Optical flash-match for this slice of TPC charge

      bool is_clear_cosmic; //!< Whether pandora marks the slice as a "clear" cosmic
      int nupdg; //!< PDG code of the neutrino candidate associated with this slice 
      float nu_score; //!< Score of how neutrino-like the slice is 
      Truth truth; //!< Truth information on the slice
      

      void setDefault();

    };
} // end namespace

#endif // SRSLICE_H
//////////////////////////////////////////////////////////////////////////////
