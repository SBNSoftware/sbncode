//  SRTrueInteraction.h
// \author  psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRTRUEINTERACTION_H
#define SRTRUEINTERACTION_H

#include "SRTrueParticle.h"

namespace caf
{
  /// The SRTrueInteraction is a representation of neutrino interaction information
  class SRTrueInteraction
  {
  public:
    SRTrueInteraction();
    ~SRTrueInteraction() {  };

    short      pdg;           ///< pdg code

    int        interaction_id; ///< ID of the truth-matched interaction
    int        inttype;       ///< Interaction type enum int_type::[...]

    bool       iscc;          ///< True if charged-current, false if not
    bool       isvtxcont;     ///< If true vertex is within TPC

    float      E;             ///< True energy [GeV]
    float      visE;          ///< True interaction deposited energy
    float      visEinslc;     ///< True deposited energy in slice [GeV]
    float      visEcosmic;    ///< True slice deposited energy from cosmics
    float      eff;           ///< Slice efficiency for this interaction
    float      pur;           ///< Slicer purity for this interaction
    float      time;          ///< interaction time.
    float      genweight;     ///< Weight, if any, assigned by the generator
    float      xsec;          ///< xsec for thrown interaction, in 1/GeV^2, as stored by the GENIE spline
    float      q2;            ///< Squared momentum transfer [GeV^2]
    float      x;             ///< Bjorken x = (k-k')^2/(2*p.q) [Dimensionless]
    float      y;             ///< Bjorken y =  (p.q)/(k.p), fractional energy loss of incoming particle [Dimensionless]
    float      w2;            ///< Final state invariant mass squared [GeV^2]

    TVector3   p;             ///< True momentum [GeV]
    TVector3   vtx;           ///< Vertex position in detector coord. [cm]

    Det_t             det;

    interaction_mode_ mode;       ///< True mode of from enum

    generator_        generator;  ///< The generator that created this neutrino interaction
    std::vector<unsigned int>   genVersion; ///< Version of the generator that created this neutrino interaction

    std::vector<SRTrueParticle> prim;       ///< Primary daughters, lepton comes first in vector.

  };

} // end namespace

#endif // SRTRUEINTERACTION_H
//////////////////////////////////////////////////////////////////////////////
