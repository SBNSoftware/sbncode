/**
 * \file Config.h
 * \brief Configuration for the LEE truth selection
 * \author A. Mastbaum <mastbaum@uchicago.edu>
 */

#ifndef __sbnanalysis_ana_LEETruthSelection_Config__
#define __sbnanalysis_ana_LEETruthSelection_Config__

#include "Util.h"
#include "MisID.h"

namespace Json {
  class Value;
}

namespace ana {
  namespace LEETruthSelection {

/** Event categories. */
enum EventType {
  k0p, k1p, kNp, kNtrk, kAny, kNSelections
};


/**
 * \class Config
 * \brief Configuration management.
 */
class Config {
public:
  /** Constructor. */
  Config() {}

  /**
   * Constructor.
   *
   * Populates all configuation parameters from the JSON settings.
   *
   * \param config The configuration as a Json object.
   */
  Config(Json::Value* config);

  /**
   * Load from JSON object.
   *
   * \param config The JSON configuration
   */
  void Initialize(Json::Value* config);

  /** Global options */
  int ntrials;  //!< Number of random trials to accept/reject events
  int dataset_id;  //!< An identifier written to the output tree
  float track_energy_distortion;  //!< Amount to smear the track energy
  float shower_energy_distortion;  //!< Amount to smear the shower energy

  /** Particle ID Matrix (binned by energy) */
  ana::LEETruthSelection::EnergyMap<ana::LEETruthSelection::PDGConfusionMatrix> pdgid_matrix;

  /** Selection types */
  bool accept_1l1p;  //!< Enable 1l1p selection
  bool accept_1l0p;  //!< Enable 1l0p selection
  bool accept_1lnp;  //!< Enable 1lNp selection
  bool accept_1lntrk;  //!< Enable 1lNtracks selection
  std::vector<EventType> selections;  //!< Enable selection types

  /** Producers for data products. */
  std::string flux_weight_producer;  //!< Flux weight producer
  std::string event_weight_producer;  //!< MC event weight weight producer
  std::string mctruth_producer;  //!< MC truth producer
  std::string mcshower_producer;  //!< MC shower producer
  std::string mctrack_producer;  //!< MC track producer
};

  }  // namespace LEETruthSelection
}  // namespace ana

#endif  // __sbnanalysis_ana_LEETruthSelection_Config__

