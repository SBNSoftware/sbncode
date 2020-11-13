//  SRTrueInteraction.h
// \author  psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRTRUEINTERACTION_H
#define SRTRUEINTERACTION_H

#include "SRTruthMatch.h"
#include "SRTrueParticle.h"

namespace caf
{
  /// The SRTrueInteraction is a representation of neutrino interaction information
  class SRTrueInteraction
  {
  public:
    SRTrueInteraction();
    ~SRTrueInteraction() {  };

    int   initpdg;         //!< Initial PDG code of probe neutrino
    int   pdg;             //!< PDG code of probe neutrino
    int   inttype;         //!< Interaction type enum int_type::[...]
    int   index;           //!< Index of the matched true neutrino interaction (-1 if not matched to neutrino)
    int   targetPDG;       //!< PDG code of struck target
    int   genie_intcode;   //!< Interaction mode (as for LArSoft MCNeutrino::Mode() )
    int   parentPDG;       //!< Parent hadron/muon PDG
    int   parentDecayMode; //!< Parent hadron/muon decay mode

    bool    isnc;              //!< same as LArSoft "ccnc" - 0=CC, 1=NC
    bool    iscc;              //!< CC (true) or NC/interference (false)
    bool    isvtxcont;         //!< If true vertex is within TPC
    bool    is_numucc_primary; //!< Whether this is the "primary" reco neutrino slice as defined by the numu CC analysis

    float      E;             ///< True energy [GeV]
    float      visE;          ///< True interaction deposited energy
    float      visEinslc;     ///< True deposited energy in slice [GeV]
    float      visEcosmic;    ///< True slice deposited energy from cosmics
    float      time;           ///< Time
    float      bjorkenX;          //!< Bjorken x = (k-k')^2/(2*p.q) [Dimensionless]
    float      inelasticityY;     //!< Inelasticity y
    float      Q2;                //!< Q squared
    float      q0;                //!< q0, struck nucleon rest frame
    float      modq;              //!< |q|, struck nucleon rest frame
    float      q0_lab;            //!< q0, lab frame
    float      modq_lab;          //!< |q|, lab frame
    float      w;                 //!< Hadronic invariant mass W
    float      t;                 //!< Kinematic t
    float      eccqe;             //!< CCQE energy
    float      baseline;          //!< Distance from decay to interaction

    SRVector3D        vtx;             //!< Vertex position in detector coord. [cm]
    SRVector3D        momentum;        //!< Neutrino three-momentum
    SRVector3D        position;        //!< Neutrino interaction position
    SRVector3D        parentDecayVtx;  //!< Parent hadron/muon decay vertex

    Det_t             det;

    interaction_mode_ mode;       ///< True mode of from enum

    generator_        generator;  ///< The generator that created this neutrino interaction
    std::vector<unsigned int>   genVersion; ///< Version of the generator that created this neutrino interaction

    std::vector<SRTrueParticle> prim;       ///< Primary daughters, lepton comes first in vector.

  };

} // end namespace

#endif // SRTRUEINTERACTION_H
//////////////////////////////////////////////////////////////////////////////
