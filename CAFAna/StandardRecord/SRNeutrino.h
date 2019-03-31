////////////////////////////////////////////////////////////////////////
// \author  Bruno Zamorano
// \date    February 2018
////////////////////////////////////////////////////////////////////////
#ifndef SRNEUTRINO_H
#define SRNEUTRINO_H

#include "TVector3.h"

namespace caf
{
  /// The SRNeutrino is a representation of neutrino interaction information
  class SRNeutrino
  {
  public:
    SRNeutrino();
    ~SRNeutrino() {  }

    bool isnc;                //!< same as LArSoft "ccnc" - 0=CC, 1=NC
    bool iscc;                //!< CC (true) or NC/interference (false)
    int pdg;                  //!< PDG code of probe neutrino
    int targetPDG;            //!< PDG code of struck target
    int genie_intcode;        //!< Interaction mode (as for LArSoft MCNeutrino::Mode() )
    double bjorkenX;          //!< Bjorken x
    double inelasticityY;     //!< Inelasticity y
    double Q2;                //!< Q squared
    double q0;                //!< q0, struck nucleon rest frame
    double modq;              //!< |q|, struck nucleon rest frame
    double q0_lab;            //!< q0, lab frame
    double modq_lab;          //!< |q|, lab frame
    double w;                 //!< Hadronic invariant mass W
    double t;                 //!< Kinematic t
    double eccqe;             //!< CCQE energy
    double energy;            //!< Neutrino energy (GeV)
    TVector3 momentum;        //!< Neutrino three-momentum
    TVector3 position;        //!< Neutrino interaction position
    int parentPDG;            //!< Parent hadron/muon PDG
    int parentDecayMode;      //!< Parent hadron/muon decay mode
    TVector3 parentDecayVtx;  //!< Parent hadron/muon decay vertex
  };

} // end namespace

#endif // SRNEUTRINO_H
//////////////////////////////////////////////////////////////////////////////
