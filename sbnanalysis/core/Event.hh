#ifndef __sbnanalysis_core_Event__
#define __sbnanalysis_core_Event__

/**
 * \file Event.hh
 *
 * The standard minimum output tree.
 *
 * This event structure is written out by every Processor subclass.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/01/25
 */

#include <map>
#include <string>
#include <vector>
#include <TTree.h>
#include <TVector3.h>

/**
 * \class Event
 * \brief The standard event data definition.
 */
class Event {
public:
  /**
   * \class Event::Metadata
   * \brief Event-level information
   */
  class Metadata {
  public:
    int run;  //!< Run ID
    int subrun;  //!< Subrun ID
    int eventID;  //!< Event ID
  };

  /**
   * \class Event::Neutrino
   * \brief Neutrino interaction information
   */
  class Neutrino {
  public:
    bool isnc;             //!< same as LArSoft "ccnc" - 0=CC, 1=NC
    bool iscc;             //!< CC (true) or NC/interference (false)
    int pdg;               //!< PDG code of probe neutrino
    int targetPDG;         //!< PDG code of struck target
    int genie_intcode;     //!< Interaction mode (as for LArSoft MCNeutrino::Mode() )
    double bjorkenX;       //!< Bjorken x
    double inelasticityY;  //!< Inelasticity y
    double Q2;             //!< Q squared
    double q0;             //!< q0, struck nucleon rest frame
    double modq;           //!< |q|, struck nucleon rest frame
    double q0_lab;         //!< q0, lab frame
    double modq_lab;       //!< |q|, lab frame
    double w;              //!< Hadronic invariant mass W
    double t;              //!< Kinematic t
    double energy;         //!< Neutrino energy (GeV)
    TVector3 momentum;     //!< Neutrino three-momentum
    
    static const int kUnfilled = -999; //!< Value all members are initialised to.
    Neutrino() : isnc(0), iscc(0),
     pdg(kUnfilled),targetPDG(kUnfilled),genie_intcode(kUnfilled),
     bjorkenX(kUnfilled),inelasticityY(kUnfilled),Q2(kUnfilled),
     q0(kUnfilled),modq(kUnfilled),q0_lab(kUnfilled),modq_lab(kUnfilled),
     w(kUnfilled),t(kUnfilled),energy(kUnfilled),
     momentum(kUnfilled,kUnfilled,kUnfilled) {};
  };

  /**
   * \class Event::FinalStateParticle
   * \brief Final state particle information
   */
  class FinalStateParticle {
  public:
    int pdg;  //!< PDG Code
    double energy;  //!< Energy
    TVector3 momentum;  //!< Three-momentum
  };

  /**
   * \class Event::Interaction
   * \brief All truth information associated with one neutrino interaction
   */
  class Interaction {
  public:
    Neutrino neutrino;  //!< The neutrino
    FinalStateParticle lepton;  //!< The primary final state lepton
    std::vector<FinalStateParticle> finalstate;  //!< The other final state particles
  };

  Metadata metadata;  //!< Event metadata
  std::vector<Interaction> interactions;  //!< All interactions

  /**
   * Event weights.
   *
   * This is a map from the weight calculator name to the list of weights
   * for all the sampled universes.
   */
  std::map<std::string, std::vector<double> > weights;
};

#endif  // __sbnanalysis_core_Event__

