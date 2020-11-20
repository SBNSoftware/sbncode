#include "CAFAna/Core/Binning.h"
#include "SBNAna/Vars/Binnings.h"


// Assume these fidutial volumes for plotting purposes, different from active volume boundaries
// x boundary reflected on the cathode atm so we'll use abs values
// {xmin = -175, xmax = -1.5}, {xmin = 1.5, xmax = 175}
std::map<std::string,double> fvnd =
{
	{"xmin", 20.},
	{"xmax", 180.},
	{"ymin", -200.},
	{"ymax", 200.},
	{"zmin", 50.},
	{"zmax", 450.}
};

// icarus assumed fidutial volume cryo 1
std::map<std::string,double> fvfd_cryo1 =
{
	{"xmin", -400.},
	{"xmax", -100.},
	{"ymin", -200.},
	{"ymax", 150.},
	{"zmin", 1000.},
	{"zmax", 1000.}
};

// icarus fidutial volume cryo 2 same as cryo 1 atm
std::map<std::string,double> fvfd_cryo2 = fvfd_cryo1;


namespace ana
{

  const Binning kNueEnergyBinning  	= Binning::Simple( 20, 0, 5);
  const Binning kNumuEnergyBinning  = Binning::Simple( 20, 0, 5);

  const Binning kLowEnergyBinning  = Binning::Custom({0.1, 0.2, 0.3, 0.4, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.});

  // Position binnings
  const Binning kPositionXNDBinning = Binning::Simple(20, fvnd["xmin"], fvnd["xmax"]);
  const Binning kPositionYNDBinning = Binning::Simple(40, fvnd["ymin"], fvnd["ymax"]);
  const Binning kPositionZNDBinning = Binning::Simple(50, fvnd["zmin"], fvnd["zmax"]);
  const Binning kPositionXFDBinning = Binning::Simple(50, fvfd_cryo1["xmin"], fvfd_cryo1["xmax"]);
  const Binning kPositionYFDBinning = Binning::Simple(35, fvfd_cryo1["ymin"], fvfd_cryo1["ymax"]);
  const Binning kPositionZFDBinning = Binning::Simple(200, fvfd_cryo1["zmin"], fvfd_cryo1["zmax"]);

}