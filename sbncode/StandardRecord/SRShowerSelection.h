////////////////////////////////////////////////////////////////////////
// \file    SRShowerSelection.h
////////////////////////////////////////////////////////////////////////
#ifndef SRSHOWERSELECTION_H
#define SRSHOWERSELECTION_H

#include "sbncode/StandardRecord/SRVector3D.h"

namespace caf
{
  /// Shower Selection metrics calculated by ShowerSelectionVals
  class SRShowerSelection
    {
    public:
      SRShowerSelection();
      ~SRShowerSelection(){  }

      // density gradient metics: split the shower up into N segments,
      // fit to the form of [0]/(X^[1]) and extract the fit
      float densityGradient;      ///< Constant in the density gradient fit
      float densityGradientPower; ///< Power in the density gradient fit

      // shower track fit metrics: fit a recob::Track to the shower stub
      // using Pandora sliding linear fit, then extract variabls
      float trackLength; ///< Lenth of fitted track [cm]
      float trackWidth;  ///< Width of fitted track (Average redidual) [cm]

      // shower residuals: DCA of each other shower in slice to the slice in question
      std::vector<float> showerResiduals; ///< Vector of residuals, size (sliceShowers-1) [cm]
    };

} // end namespace

#endif // SRSHOWERSELECTION_H
//////////////////////////////////////////////////////////////////////////////
