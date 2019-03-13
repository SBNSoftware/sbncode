////////////////////////////////////////////////////////////////////////
// \author  Bruno Zamorano
// \date    December 2018
////////////////////////////////////////////////////////////////////////
#ifndef SRHEADERBRANCH_H
#define SRHEADERBRANCH_H

#include <vector>

namespace caf
{
  /// \brief Contains information about the run as well as file metadata
  class SRHeaderBranch
    {
    public:
      SRHeaderBranch();
      ~SRHeaderBranch();

      int          run;           ///< Run number
      int       subrun;           ///< Subrun number
      int      eventID;           ///< Event ID number

      void setDefault();

    };
  
} // end namespace

#endif // SRHEADERBRANCH_H
//////////////////////////////////////////////////////////////////////////////
