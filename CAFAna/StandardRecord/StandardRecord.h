#ifndef STANDARDRECORD_H
#define STANDARDRECORD_H

#include "StandardRecord/SRMetadata.h"
#include "StandardRecord/SRInteraction.h"
#include "StandardRecord/SRRecoInteraction.h"

/// Common Analysis Files
namespace caf
{
  
  /// \brief   The StandardRecord is the primary top-level object in the 
  ///          Common Analysis File trees.   
  
  // Based on sbnanalysis Event
  class StandardRecord  
  {
  public:
    StandardRecord();
    ~StandardRecord();

    SRMetadata metadata;  //!< Event metadata
    std::vector<SRInteraction> truth; //!< All truth interactions
    std::vector<SRRecoInteraction> reco; //!< Reconstructed interactions

    //    Experiment experiment;  //!< Experiment identifier
  };
  
} // end namespace

#endif // STANDARDRECORD_H
//////////////////////////////////////////////////////////////////////////////
