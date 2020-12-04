#include "SBNAna/Cuts/VolumeDefinitions.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include <cmath>

bool PtInVol(const caf::SRVector3DProxy& pt, const FidVol& vol)
{
  return (vol.xmin < pt.x && pt.x < vol.xmax &&
          vol.ymin < pt.y && pt.y < vol.ymax &&
          vol.zmin < pt.z && pt.z < vol.zmax);
}


bool PtInVolAbsX(const caf::SRVector3DProxy& pt, const FidVol& vol)
{
  return (vol.xmin < fabs(pt.x) && fabs(pt.x) < vol.xmax &&
          vol.ymin <      pt.y  &&      pt.y  < vol.ymax &&
          vol.zmin <      pt.z  &&      pt.z  < vol.zmax);
}
