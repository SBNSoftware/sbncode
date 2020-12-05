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
  
  const std::vector<double> kLowEnergyGeVEdges = {.1, .2, .3, .4, .5, .75,
						    1., 1.25, 1.5, 1.75,
						    2., 2.25, 2.5, 2.75,
						    3., 3.25, 3.5, 3.75,
						    4., 4.25, 4.5, 4.75, 5.};

  // Take note on the dimensitions of the variables instored in the cafs
  const Binning kNueEnergyBinning   = Binning::Simple( 20, 0., 5000.); // MeV
  const Binning kNumuEnergyBinning  = Binning::Simple( 20, 0., 5000.); // MeV

  const Binning kLowEnergyBinning     = Binning::Custom( kLowEnergyEdges ); // MeV
  const Binning kLowEnergyGeVBinning  = Binning::Custom( kLowEnergyGeVEdges ); // GeV

  // Position binnings
  const Binning kPositionXNDBinning = Binning::Simple(40, avnd.xmin, avnd.xmax);
  const Binning kPositionYNDBinning = Binning::Simple(40, avnd.ymin, avnd.ymax);
  const Binning kPositionZNDBinning = Binning::Simple(50, avnd.zmin, avnd.zmax);
  const Binning kPositionXFDBinning = Binning::Simple(47, avfd_cryo1.xmin, avfd_cryo1.xmax);
  const Binning kPositionYFDBinning = Binning::Simple(35, avfd_cryo1.ymin, avfd_cryo1.ymax);
  const Binning kPositionZFDBinning = Binning::Simple(90, avfd_cryo1.zmin, avfd_cryo1.zmax);

  const Binning kCRTXFDBinning = Binning::Simple(60, crtfd.xmin, crtfd.xmax);
  const Binning kCRTYFDBinning = Binning::Simple(100, crtfd.ymin, crtfd.ymax);
  const Binning kCRTZFDBinning = Binning::Simple(130, crtfd.zmin, crtfd.zmax);

  const Binning kBDTBinning       = Binning::Simple(50,0, 1.0);

}
