////////////////////////////////////////////////////////////////////////
// \file    SRSlice.h
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRSLICE_H
#define SRSLICE_H

#include "SRTrueInteraction.h"
#include "SRTruthMatch.h"
#include "SRVector3D.h"


namespace caf
{
  /// A matching of TPC slice charge to Optical flash light
  class FlashMatch {
  public:
    bool present;
    float score;
    float time;
    float pe;
  };

  /// Matching between this SRSlice and the corresponding SRTrueInteraction
  class TruthMatch {
  public:
    float      visEinslc;         ///< True deposited energy in slice [GeV]
    float      visEcosmic;        ///< True slice deposited energy from cosmics
    float      eff;               ///< Slice efficiency for this interaction
    float      pur;               ///< Slicer purity for this interaction
    int        index;             ///< Index of the matched true neutrino interaction 
                                  ///< (-1 if not matched to neutrino)
    bool       is_numucc_primary; ///< Whether this is the "primary" reco neutrino slice 
                                  ///< as defined by the numu CC analysis
  };

  /// An  SRSlice contains overarching information for a slice.
  class SRSlice
    {
    public:

      SRSlice();
      virtual ~SRSlice();

      unsigned producer;    ///< Index of the producer that produced this object. 
                            ///< In ICARUS, this is the same as the cryostat.
      float    charge;      ///< Calorimetric energy
      SRVector3D vertex;    ///< Candidate neutrino vertex in local detector coordinates [cm]

      SRTrueInteraction truth; //!< Truth information on the slice      
      SRTruthMatch tmatch; //!< Matching information between truth and reco objects

      SRFlashMatch fmatch_a; //!< Optical flash-match for this slice of TPC charge
      SRFlashMatch fmatch_b; //!< Optical flash-match for this slice of TPC charge

      FlashMatch fmatch; //!< Optical flash-match for this slice of TPC charge

      bool is_clear_cosmic; //!< Whether pandora marks the slice as a "clear" cosmic
      int nu_pdg; //!< PDG assigned to the PFParticle Neutrino
      float nu_score; //!< Score of how neutrino-like the slice is 
      std::vector<size_t> primary; //!< ID's of primary tracks and showers in slice
      int                 self;    //!< ID of the particle representing this slice

      void setDefault();

    };
} // end namespace

#endif // SRSLICE_H
//////////////////////////////////////////////////////////////////////////////
