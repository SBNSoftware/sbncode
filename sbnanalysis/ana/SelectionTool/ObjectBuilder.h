#ifndef __sbnanalysis_ana_SelectionTool_ObjectBuilder__
#define __sbnanalysis_ana_SelectionTool_ObjectBuilder__

/**
 * \file ObjectBuilder.h
 *
 * A processor to build the objects for the selection tool.
 *
 * We define the methods called for initialization, finalization, 
 * and event-by-event processing.
 *
 * Author: R. Jones <rjones@hep.ph.liv.ac.uk>
 */

// Includes
#include <iostream>
#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"

// Forward declarations go here

/** All analysis code is defined in namespace "ana" */
namespace ana {

  /** Code specific to the SelectionTool. */
  namespace SelectionTool {

    /**
     * \class ObjectBuilder
     * \brief Building selection tool objects from gallery/LArSoft reco events
     */
    class ObjectBuilder {
      public:

        /** Constructor. */
        ObjectBuilder();

        /**
         * Initialization.
         *
         * Here we load configuration parameters, set up histograms for output, and
         * add our own branches to the output tree.
         *
         * \param config A configuration, as a JSON object
         */
        //void Initialize(Json::Value* config);
        
        /**
         * Initialization. (ONCE FHICL HAS BEEN MERGED)
         *
         * Here we load configuration parameters, set up histograms for output, and
         * add our own branches to the output tree.
         *
         * \param config A configuration, as a FHICL object, same as in LArSoft
         */
        void Initialize(fhicl::ParameterSet const &p);

        /** Finalize and write objects to the output file. */
        void Finalize();

        /**
         * Process one event.
         *
         * \param ev A single event, as a gallery::Event
         * \param reco Reconstructed interactions
         * \return True to keep event
         */
        bool ProcessEvent(const gallery::Event& ev, std::vector<Event::RecoInteraction>& reco);


      protected:

        // Member variables from the external fhicl file
        // Handle labels                                                                
        art::InputTag fGeneratorLabel;                                                  
        art::InputTag fGeantLabel;                                                      
        art::InputTag fPandoraLabel;                                                    
        art::InputTag fRecoTrackLabel;                                                 
        art::InputTag fRecoShowerLabel;                                                
        art::InputTag fRecoTrackCalorimetryLabel;                                     
        art::InputTag fRecoTrackParticleidLabel;                                      
        art::InputTag fHitLabel; 

    }; // ObjectBuilder class
  }  // namespace SelectionTool
}  // namespace ana

#endif  // __sbnanalysis_ana_SelectionTool_ObjectBuilder__

