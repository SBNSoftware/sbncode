#ifndef STANDARDRECORD_H
#define STANDARDRECORD_H

#include "StandardRecord/SRTruthBranch.h"

/// Common Analysis Files
namespace caf
{
  
  /// \brief   The StandardRecord is the primary top-level object in the 
  ///          Common Analysis File trees.   
  
  class StandardRecord  
  {
    
  public:
    StandardRecord();
    ~StandardRecord();

    SRTruthBranch    truth;     ///< Truth branch for MC: energy, flavor, etc.
  };
  
} // end namespace

#endif // STANDARDRECORD_H
//////////////////////////////////////////////////////////////////////////////
