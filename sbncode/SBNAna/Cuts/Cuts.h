#pragma once

// Definition of the generic Cut object
#include "CAFAna/Core/Cut.h"

// =================================================================== //
// These values have been used for the nue event selection as of       // 
// April 29 2020. The active and fiducial volumes are the same so far  //
// =================================================================== //

// x boundary reflected on the cathode atm so we'll use abs values
// {xmin = -175, xmax = -1.5}, {xmin = 1.5, xmax = 175}
std::map<std::string,double> avnd =
{
	{"xmin", 1.5},
	{"xmax", 175.},
	{"ymin", -175.},
	{"ymax", 175.},
	{"zmin", 30.},
	{"zmax", 450.}
};

// icarus active volume cryo 1
std::map<std::string,double> avfd_cryo1 =
{
	{"xmin", -368.49},
	{"xmax", -71.94},
	{"ymin", -181.86},
	{"ymax", 134.96},
	{"zmin", -894.95},
	{"zmax", 894.95}
};

// icarus active volume cryo 2 same as cryo 1 atm
std::map<std::string,double> avfd_cryo2 = avfd_cryo1;

namespace ana{

  extern const Cut kActiveVolumeND;
  extern const Cut kActiveVolumeFDCryo1;
  extern const Cut kActiveVolumeFDCryo2;

} // namespace
