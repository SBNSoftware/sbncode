#include "CAFAna/Core/Binning.h"
#include "SBNAna/Vars/Binnings.h"

#include "SBNAna/Cuts/VolumeDefinitions.h"

namespace ana
{

  const std::vector<double> kLowEnergyEdges = {100., 200., 300., 400., 500., 750.,
					       1000., 1250., 1500., 1750.,
					       2000., 2250., 2500., 2750.,
					       3000., 3250., 3500., 3750.,
					       4000., 4250., 4500., 4750., 5000.}; 
  
  const Binning kNueEnergyBinning   = Binning::Simple( 20, 0., 5000.); // MeV
  const Binning kNumuEnergyBinning  = Binning::Simple( 20, 0., 5000.); // MeV

  const Binning kLowEnergyBinning  = Binning::Custom( kLowEnergyEdges ); // MeV

  // Position binnings
  const Binning kPositionXNDBinning = Binning::Simple(40, avnd.xmin, avnd.xmax);
  const Binning kPositionYNDBinning = Binning::Simple(40, avnd.ymin, avnd.ymax);
  const Binning kPositionZNDBinning = Binning::Simple(50, avnd.zmin, avnd.zmax);
  const Binning kPositionXFDBinning = Binning::Simple(47, avfd_cryo1.xmin, avfd_cryo1.xmax);
  const Binning kPositionYFDBinning = Binning::Simple(35, avfd_cryo1.ymin, avfd_cryo1.ymax);
  const Binning kPositionZFDBinning = Binning::Simple(90, avfd_cryo1.zmin, avfd_cryo1.zmax);

}
