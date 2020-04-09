////////////////////////////////////////////////////////////////////////
// \file    SRSlice.h
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICE_H
#define SRSLICE_H

#include "SRTrueInteraction.h"
#include "SRVector3D.h"


namespace caf
{
  /// An  SRSlice contains overarching information for a slice.
  class SRSlice
    {
    public:

      /// A matching of TPC slice charge to Optical flash light
      struct FlashMatch {
        bool present;
        float score;
        float time;
        float pe;
      };

      /// Matching between this SRSlice and the corresponding SRTrueInteraction
      struct TruthMatch {
	float      visEinslc;     ///< True deposited energy in slice [GeV]
	float      visEcosmic;    ///< True slice deposited energy from cosmics
	float      eff;           ///< Slice efficiency for this interaction
	float      pur;           ///< Slicer purity for this interaction
        int        index;         ///< Index of the matched true neutrino interaction (-1 if not matched to neutrino)
        bool       is_numucc_primary; ///< Whether this is the "primary" reco neutrino slice as defined by the numu CC analysis
      };

      SRSlice();
      virtual ~SRSlice();


      int      id;          ///< number of hits
      float    charge;      ///< Calorimetric energy
      SRVector3D vertex;      ///< Candidate neutrino vertex in local detector coordinates [cm]

      SRTrueInteraction truth; //!< Truth information on the slice      
      TruthMatch tmatch; //!< Matching information between truth and reco objects

      void setDefault();

      // To Do: Make SRFlashMatch and add here
      FlashMatch fmatch; //!< Optical flash-match for this slice of TPC charge

      bool is_clear_cosmic; //!< Whether pandora marks the slice as a "clear" cosmic
      float nu_score; //!< Score of how neutrino-like the slice is 
      std::vector<size_t> primary; //!< ID's of primary tracks and showers in slice
      int                 self;    //!< ID of the particle representing this slice



    };
} // end namespace

#endif // SRSLICE_H
//////////////////////////////////////////////////////////////////////////////
