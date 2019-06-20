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
#include "Experiment.hh"

namespace event {

static const int kUnfilled = -99999;  //!< Value for unfilled variables


/**
 * \class Metadata
 * \brief Event-level information
 */
class Metadata {
public:
  /** Constructor. */
  Metadata() : run(kUnfilled), subrun(kUnfilled), eventID(kUnfilled),
               mcIndex(kUnfilled) {}

  /** Reset members to defaults. */
  void Init() {
    run = kUnfilled;
    subrun = kUnfilled;
    eventID = kUnfilled;
    mcIndex = kUnfilled;
  }

  int run;      //!< Run ID
  int subrun;   //!< Subrun ID
  int eventID;  //!< Event ID
  int mcIndex;  //!< MC event index
};


/**
 * \class Neutrino
 * \brief Neutrino interaction information
 */
class Neutrino {
public:
  /** Constructor. */
  Neutrino()
    : isnc(false), iscc(false), pdg(0), initpdg(0), targetPDG(0),
      genie_intcode(0), bjorkenX(kUnfilled), inelasticityY(kUnfilled),
      Q2(kUnfilled), q0(kUnfilled),
      modq(kUnfilled), q0_lab(kUnfilled), modq_lab(kUnfilled),
      w(kUnfilled), t(kUnfilled), energy(kUnfilled),
      momentum(kUnfilled, kUnfilled, kUnfilled), parentPDG(0),
      parentDecayMode(0), parentDecayVtx(kUnfilled, kUnfilled, kUnfilled),
      baseline(kUnfilled) {}

  bool isnc;                //!< same as LArSoft "ccnc" - 0=CC, 1=NC
  bool iscc;                //!< CC (true) or NC/interference (false)
  int initpdg;              //!< Initial PDG code of probe neutrino
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
  double baseline;          //!< Distance from decay to interaction
};


/**
 * \class FinalStateParticle
 * \brief Final state particle information
 */
class FinalStateParticle {
public:
  /** Constructor. */
  FinalStateParticle()
    : pdg(kUnfilled), energy(kUnfilled),
      momentum(kUnfilled, kUnfilled, kUnfilled) {}

  /** GENIE status codes, see Framework/GHEP/GHepStatus.h for details */
  enum GenieStatus {
    kIStUndefined                = -1, 
    kIStInitialState             =  0,
    kIStStableFinalState         =  1,
    kIStIntermediateState        =  2,
    kIStDecayedState             =  3,
    kIStCorrelatedNucleon        = 10,
    kIStNucleonTarget            = 11,
    kIStDISPreFragmHadronicState = 12,
    kIStPreDecayResonantState    = 13,
    kIStHadronInTheNucleus       = 14,
    kIStFinalStateNuclearRemnant = 15,
    kIStNucleonClusterTarget     = 16
  };

  int pdg;  //!< PDG Code
  double energy;  //!< Energy
  TVector3 momentum;  //!< Three-momentum
  TVector3 start;  //!< Start position in detector coords [cm]
  TVector3 end;  //!< End position in detector coords (may be outside of AV) [cm]

  /**
   * Length of particle contained in the detector active volume [cm]
   * Set to -1 if no Detector is defined through the ProviderManager
   * Set to 0 for a particle that originates outside the detector active volume
   */
  double contained_length;

  double length; //!< Total length of the energy depositions [cm]
  bool is_primary; //!< Whether the process producing the particle was "primary"
  int status_code; //!< Status code returned by GENIE (see GenieStatus enum)
};


/**
 * \class Interaction
 * \brief All truth information associated with one neutrino interaction
 */
class Interaction {
public:
  Neutrino neutrino;  //!< The neutrino
  FinalStateParticle lepton;  //!< The primary final state lepton
  std::vector<FinalStateParticle> finalstate; //!< Other final state particles
  size_t nfinalstate;  //!< Size of finalstate

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
    RecoInteraction()
        : truth_index(-1), reco_energy(kUnfilled), weight(1) {}

    /** Fill in truth information -- other fields set as in default */
    explicit RecoInteraction(const Interaction &t, int index)
        : truth(t), truth_index(index), reco_energy(kUnfilled), weight(1) {}

    Interaction truth;  //!< Truth level information about interaction

    /**
     * Index into the vector of truth interaction objects in the Event
     * (same as the index into MCTruth objects). Equal to -1 if there is
     * no corresponding truth interaction.
     */
    int truth_index;

    double reco_energy;  //!< Reconstructed neutrino energy [GeV]
    double weight;  //!< Selection-defined event weight
};


/**
 * \class RecoEvent
 * \brief The reconstructed event data definition.
 *
 * An Event contains a single reconstructed neutrino with the corresponding
 * truth information for the best-matched true neutrino, if available.
 */
class RecoEvent {
public:
  Metadata metadata;  //!< Event metadata
  RecoInteraction reco; //!< Reconstructed interactions
  Experiment experiment;  //!< Experiment identifier
};


/**
 * \class Event
 * \brief The standard event data definition.
 *
 * An Event contains all Monte Carlo events, with lists of all true
 * and reconstructed neutrinos.
 */
class Event {
public:
  Metadata metadata;  //!< Event metadata
  size_t ntruth;  //!< Size of truth
  std::vector<Interaction> truth; //!< All truth interactions
  size_t nreco;  //!< Size of reco
  std::vector<RecoInteraction> reco; //!< Reconstructed interactions
  Experiment experiment;  //!< Experiment identifier
};

}  // namespace event

#endif  // __sbnanalysis_core_Event__

