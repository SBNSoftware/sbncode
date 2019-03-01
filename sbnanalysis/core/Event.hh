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
    /** Constructor. */
    Metadata() : run(kUnfilled), subrun(kUnfilled), eventID(kUnfilled) {}

    /** Reset members to defaults. */
    void Init() {
      run = kUnfilled;
      subrun = kUnfilled;
      eventID = kUnfilled;
    }

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
    /** Constructor. */
    Neutrino()
      : isnc(false), iscc(false), pdg(0), targetPDG(0), genie_intcode(0),
        bjorkenX(kUnfilled), inelasticityY(kUnfilled), Q2(kUnfilled),
        q0(kUnfilled), modq(kUnfilled), q0_lab(kUnfilled), modq_lab(kUnfilled),
        w(kUnfilled), t(kUnfilled), energy(kUnfilled),
        momentum(kUnfilled, kUnfilled, kUnfilled) {}

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
    double eccqe;          //!< CCQE energy
    double energy;         //!< Neutrino energy (GeV)
    TVector3 momentum;     //!< Neutrino three-momentum
    TVector3 position;     //!< Neutrino interaction position
  };

  /**
   * \class Event::FinalStateParticle
   * \brief Final state particle information
   */
  class FinalStateParticle {
  public:
    /** Constructor. */
    FinalStateParticle()
      : pdg(kUnfilled), energy(kUnfilled),
        momentum(kUnfilled, kUnfilled, kUnfilled) {}

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

    /** The other final state particles. */
    std::vector<FinalStateParticle> finalstate; //!< Final state particles

    /**
     * Event weights.
     *
     * This is a map from the weight calculator name to the list of weights
     * for all the sampled universes.
     */
    std::map<std::string, std::vector<double> > weights;
  };

  /**
   * \class RecoInteraction
   * \brief Contains truth level information and additional fields for
   * selection-defined reconstruction information
   */
  class RecoInteraction {
    public:
      /** Default Constructor */
      RecoInteraction(): 
        truth_index(-1), 
        reco_energy(kUnfilled),
        weight(1.) {}

      /** Fill in truth information -- other fields set as in default */
      explicit RecoInteraction(const Interaction &t, int index): 
        truth(t), 
        truth_index(index),
        reco_energy(kUnfilled),
        weight(1.) {}

      Interaction truth; //!< Contains truth level information about interaction

      /**
       * Index into the vector of truth interaction objects in the Event
       * (same as the index into MCTruth objects). Equal to -1 if there is
       * no corresponding truth interaction.
       */
      int truth_index;

      /**
       * Selection defined reconstructed energy of neutrino. Units in GeV to keep
       * consistent w/ Interaction class. */
      double reco_energy;

      /**
       * Selection defined weight of reconstructed interaction to be used by downstream
       * analyis. */
      double weight;  
  };

  Metadata metadata;  //!< Event metadata
  std::vector<Interaction> truth; //!< All truth interactions
  std::vector<RecoInteraction> reco; //!< Reconstructed interactions

  static const int kUnfilled = -99999;  //!< Value for unfilled variables
};

#endif  // __sbnanalysis_core_Event__

