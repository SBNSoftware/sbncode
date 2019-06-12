#ifndef __sbnanalysis_ana_SelectionAnalysis_AnalysisMain__
#define __sbnanalysis_ana_SelectionAnalysis_AnalysisMain__

/**
 * \file AnalysisMain.h
 *
 * An example event selection analyser.
 *
 * This is an implementation of the core::SelectionBase class. We define
 * the methods called for initialization, finalization, and event-by-event
 * processing.
 *
 * Author: R. Jones <rjones@hep.ph.liv.ac.uk>
 */

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/Experiment.hh"

/** All analysis code is defined in namespace "ana" */
namespace ana {

  /** Code specific to the SelectionAnalysis. */
  namespace SelectionAnalysis {

    /**
     * @class AnalysisMain
     * @brief An example analysis with a selection
     *
     */
    class AnalysisMain : public core::SelectionBase {
    
      public:
        /** Constructor. */
        AnalysisMain();

        /**
         * Initialization.
         *
         * Here we load configuration parameters, set up histograms for output, and
         * add our own branches to the output tree.
         *
         * \param config A configuration, as FHiCL ParameterSet object
         */
        void Initialize(fhicl::ParameterSet* config=NULL);

        /** Finalize and write objects to the output file. */
        void Finalize();

        /**
         * Process one event.
         *
         * \param ev A single event, as a gallery::Event
         * \param reco Reconstructed interactions
         * \return True to keep event
         */
        bool ProcessEvent(const gallery::Event& ev,
                          const std::vector<Event::Interaction>& truth,
                          std::vector<Event::RecoInteraction>& reco);

      protected:
        unsigned fEventCounter;      //!< Count processed events
        unsigned fTruthEntryCounter; //!< Count processed truth interactions
        unsigned fRecoEntryCounter;  //!< Count processed reco interactions

    }; // AnalysisMain class
  }  // namespace SelectionAnalysis
}  // namespace ana

#endif  // __sbnanalysis_ana_SelectionAnalysis_AnalysisMain__

