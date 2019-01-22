#ifndef __sbnanalysis_util_SelectionTool_ObjectBuilder__
#define __sbnanalysis_util_SelectionTool_ObjectBuilder__

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

/** All analysis code is defined in namespace "util" */
namespace util {

  /** Code specific to the SelectionTool. */
  namespace SelectionTool {

    /**
     * \class ObjectBuilder
     * \brief Building selection tool objects from gallery/LArSoft reco events
     */
    class ObjectBuilder : public core::SelectionBase {
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
        void Initialize(fhicl::ParameterSet *p);

        /** Finalize and write objects to the output file. */
        void Finalize();

        /**
         * Process one event.
         *
         * \param ev A single event, as a gallery::Event
         * \return True to keep event
         */
        bool ProcessEvent(const gallery::Event& ev);


      protected:

        // Member variables from the external fhicl file
        // Handle labels                                                                
        std::string fGeneratorLabel;                                                  
        std::string fGeantLabel;                                                      
        std::string fPandoraLabel;                                                    
        std::string fRecoTrackLabel;                                                 
        std::string fRecoShowerLabel;                                                
        std::string fRecoTrackCalorimetryLabel;                                     
        std::string fRecoTrackParticleidLabel;                                      
        std::string fHitLabel; 

    }; // ObjectBuilder class
  }  // namespace SelectionTool
}  // namespace util

#endif  // __sbnanalysis_util_SelectionTool_ObjectBuilder__

