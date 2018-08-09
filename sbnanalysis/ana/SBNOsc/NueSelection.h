#ifndef __sbnanalysis_ana_SBNOsc_NueSelection__
#define __sbnanalysis_ana_SBNOsc_NueSelection__

/**
 * \file NueSelection.h
 *
 * SBN nue selection.
 *
 * Author: 
 */

#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"

class TH2D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NueSelection
 * \brief Electron neutrino event selection
 */
class NueSelection : public core::SelectionBase {
public:
  /** Constructor. */
  NueSelection();

  /**
   * Initialization.
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
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, std::vector<Event::Interaction>& reco);

protected:
  unsigned fEventCounter;  //!< Count processed events
  unsigned fNuCount;  //!< Count selected events

  /** Configuration parameters */
  art::InputTag fTruthTag;  //!< art tag for MCTruth information
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NueSelection__

