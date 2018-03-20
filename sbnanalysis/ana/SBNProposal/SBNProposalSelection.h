#ifndef __sbnanalysis_ana_SBNProposalAnalysis_SBNProposalSelection__
#define __sbnanalysis_ana_SBNProposalAnalysis_SBNProposalSelection__

/**
 * \file SBNProposalSelection.h
 *
 * This module is intended to mimic the SBN Proposal selection
 *
 * This is intended to perform:
 *  - nue selection
 *      > Based on smeared energy 
 *  - numu selection
 *      > Based on smeared muon energy
 *
 * Author: J. Zennamo <jaz8600@fnal.gov>
 */

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"

// Forward declarations
class TH2D;

/** All analysis code is defined in namespace "ana" */
namespace ana {

  /** Code specific to the SBNProposalAnalysis. */
  namespace SBNProposalAnalysis {

/**
 * \class SBNProposalSelection
 * \brief An example selection analysis
 *
 * This selection analysis doesn't actually select events, it just
 * demonstrates the framework!
 */
class SBNProposalSelection : public core::SelectionBase {
public:
  /** Constructor. */
  SBNProposalSelection();

  /**
   * Initialization.
   *
   * Here we load configuration parameters, set up histograms for output, and
   * add our own branches to the output tree.
   *
   * \param config A configuration, as a JSON object
   */
  void Initialize(Json::Value* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \return True to keep event
   */
  bool ProcessEvent(gallery::Event& ev);

protected:
  unsigned fEventCounter;  //!< Count processed events

  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
  int fMyParam;  //!< A parameter from the configuration file

  /** Custom data branches */
  int fNuCount;  //!< Number of neutrino interactions in the event
  int fMyVar;  //!< Another variable of interest

  /** Histograms */
  TH2D* fNuVertexXZHist;  //!< Neutrino vertex XZ projection
};

  }  // namespace SBNProposalAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNProposalAnalysis_SBNProposalSelection__

