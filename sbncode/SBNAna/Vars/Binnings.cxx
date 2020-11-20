#include "CAFAna/Core/Binning.h"
#include "SBNAna/Vars/Binnings.h"

#include "SBNAna/Cuts/VolumeDefinitions.h"

namespace ana
{

  const Binning kNueEnergyBinning  	= Binning::Simple( 20, 0, 5);
  const Binning kNumuEnergyBinning  = Binning::Simple( 20, 0, 5);

  const Binning kLowEnergyBinning  = Binning::Custom({0.1, 0.2, 0.3, 0.4, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.});

  // Position binnings
  const Binning kPositionXNDBinning = Binning::Simple(40, avnd.xmin, avnd.xmax);
  const Binning kPositionYNDBinning = Binning::Simple(40, avnd.ymin, avnd.ymax);
  const Binning kPositionZNDBinning = Binning::Simple(50, avnd.zmin, avnd.zmax);
  const Binning kPositionXFDBinning = Binning::Simple(50, avfd_cryo1.xmin, avfd_cryo1.xmax);
  const Binning kPositionYFDBinning = Binning::Simple(35, avfd_cryo1.ymin, avfd_cryo1.ymax);
  const Binning kPositionZFDBinning = Binning::Simple(200, avfd_cryo1.zmin, avfd_cryo1.zmax);

}
