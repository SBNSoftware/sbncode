#pragma once

#include "StandardRecord/Proxy/FwdDeclare.h"

// =================================================================== //
// These values have been used for the nue event selection as of       //
// April 29 2020. The active and fiducial volumes are the same so far  //
// =================================================================== //

struct FidVol
{
  float xmin, xmax, ymin, ymax, zmin, zmax;
};

// Only for use with FD volume
bool PtInVol(const caf::SRVector3DProxy& pt, const FidVol& vol);

// Compare (|x|, y, z) to the volume definition. Must use this for ND volume.
bool PtInVolAbsX(const caf::SRVector3DProxy& pt, const FidVol& vol);

// x boundary reflected on the cathode atm so we'll use abs values
// {xmin = -175, xmax = -1.5}, {xmin = 1.5, xmax = 175}
const FidVol fvndAbs{  +1.5, +175.,  // x
                      -175., +175.,  // y
                       +30., +450.}; // z

// Sometimes (e.g. plotting) we want the full volume 
const FidVol fvnd{ -175., +175.,  // x
                   -175., +175.,  // y
                    +30., +450.}; // z

const FidVol avnd{ -200., +200.,  // x
                   -200., +200.,  // y
                     +0., +500.}; // z

// icarus fiducial volume cryo used to be the same as the actual active volume (not the values below)
// slightly change these values. still not uptimized.
const FidVol fvfd_cryo1{ -368.49,  -71.94,  // x
                         -181.86, +134.96,  // y
                         -894.95, +894.85}; // z

// icarus cheat active volume cryo 1 for plotting
const FidVol avfd_cryo1{ -400.,  -50.,  // x
                         -200., +150.,  // y
                         -1000., +1000.}; // z

// icarus active volume cryo 2 same as cryo 1 atm
const FidVol fvfd_cryo2 = fvfd_cryo1;
const FidVol avfd_cryo2 = avfd_cryo1;

// not a fiducial volume but for plotting CRT hit positions
const FidVol crtfd{ -600,   0.,  // x
                    -300., +700.,  // y
                    -1000, +1600}; // z


// Once we have C++20 we will be able to rewrite these initializers slightly
// more explicitly like this
//
// const FidVol avnd{.xmin = +1.5, .xmax = +175, ...
