#pragma once

// =================================================================== //
// These values have been used for the nue event selection as of       //
// April 29 2020. The active and fiducial volumes are the same so far  //
// =================================================================== //

struct FidVol
{
  float xmin, xmax, ymin, ymax, zmin, zmax;
};

// x boundary reflected on the cathode atm so we'll use abs values
// {xmin = -175, xmax = -1.5}, {xmin = 1.5, xmax = 175}
const FidVol avnd{  +1.5, +175.,  // x
                  -175.,  +175.,  // y
                   +30.,  +450.}; // z

// icarus active volume cryo 1
const FidVol avfd_cryo1{ -368.49,  -71.94,  // x
                         -181.86, +134.96,  // y
                         -894.95, +894.95}; // z

// icarus active volume cryo 2 same as cryo 1 atm
const FidVol avfd_cryo2 = avfd_cryo1;


// Once we have C++20 we will be able to rewrite these initializers slightly
// more explicitly like this
//
// const FidVol avnd{.xmin = +1.5, .xmax = +175, ...
