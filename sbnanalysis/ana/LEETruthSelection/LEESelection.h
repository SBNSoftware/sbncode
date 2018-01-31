/**
 * \file LEESelection.h
 * \brief A truth-based LEE event selection
 * \author A. Mastbaum, G. Putnam
 */

#ifndef __sbnanalysis_ana_LEETruthSelection_LEESelection__
#define __sbnanalysis_ana_LEETruthSelection_LEESelection__

#include <string>
#include <vector>
#include "gallery/Event.h"
#include "core/SelectionBase.hh"

#include "Config.h"
#include "Util.h"

namespace ana {
  namespace LEETruthSelection {

/**
 * \class LEESelection
 * \brief Truth-based low-energy excess event selection
 *
 * A truth-based selection for low-energy neutrino interactions, for
 * Monte Carlo low-energy excess sensitivity studies.
 *
 * Configurable for different selection modes, with one lepton (e or
 * mu) plus no protons (1l0p), exactly one proton (1l1p), one or more
 * protons (1lnp) or any number of tracks (1lntrk).
 */
class LEESelection : public core::SelectionBase {
public:
  /** Constructor. */
  LEESelection() : SelectionBase() {}

  /**
   * Initialization.
   *
   * Load configuration and add selection output to the event tree.
   *
   * \param config A configuration, as a JSON object
   */
  void Initialize(Json::Value* config=NULL);

  /** Finalize and print summary information. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   */
  void ProcessEvent(gallery::Event& ev);

  float nextTrackEnergyDistortion(float this_energy);
  float nextShowerEnergyDistortion(float this_energy);
  int nextParticleID(float energy, int true_pdgid);

  /**
   * \struct OutputData
   * \brief A container for selection output
   */
  struct OutputData {
    /** Event information. */
    int np;  //!< Number of protons
    int ntrk;  //!< Number of tracks
    double bnbweight;  //!< BNB correction weight
    int dataset;  //!< Dataset ID

    /** Lepton */
    int lpid;  //!< Lepton PID PDG
    int lpdg;  //!< Lepton true PDG
    bool lexit;  //!< Lepton exiting?
    double levis;  //!< Lepton visible energy
    double llen;  //!< Lepton track length
    TLorentzVector lmomentum;  //!< Lepton four-momentum
    double eccqe;  //!< CCQE "reconstructed" neutrino energy

    /** Tracks */
    std::vector<int> track_pdg;  //!< Track PID PDGs
    std::vector<int> track_pdg_true;  //!< Track true PDGs
    std::vector<double> track_evis;  //!< Track energies
    std::vector<TLorentzVector> track_momentum;  //!< Track four-momentum
  };

protected:
  Config fConfig;  //!< Configuration manager
  OutputData fOutputData;  //!< Selection output container
  //EnergyMap<PDGConfusionMatrix> fMisIDMap;  //!< Particle mis-ID mapping
};

  }  // namespace LEETruthSelection
}  // namespace ana

#endif  // __sbnanalysis_ana_LEETruthSelection_LEESelection__

