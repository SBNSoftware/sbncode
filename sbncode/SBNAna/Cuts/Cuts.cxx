#include "SBNAna/Cuts/Cuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include <map>
#include <string>

// =================================================================== //
// These values have been used for the nue event selection as of       // 
// April 29 2020. The active and fiducial volumes are the same so far  //
// =================================================================== //

// x boundary reflected on the cathode atm so we'll use abs values
// {xmin = -175, xmax = -1.5}, {xmin = 1.5, xmax = 175}
std::map<std::string,double> avnd =
{
	{"xmin",1.5},
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

  const Cut kActiveVolumeND(
  	[](const caf::SRProxy* sr){
  		bool x = (avnd["xmin"] < abs(sr->slc[0].vertex.x)) && (abs(sr->slc[0].vertex.x) < avnd["xmax"]);
  		bool y = (avnd["ymin"] < sr->slc[0].vertex.y) && (sr->slc[0].vertex.y < avnd["ymax"]);
  		bool z = (avnd["zmin"] < sr->slc[0].vertex.z) && (sr->slc[0].vertex.z < avnd["zmax"]);
  		return(x && y && z);
  	}
  	);

  const Cut kActiveVolumeFDCryo1(
  	[](const caf::SRProxy* sr){
  		bool x = (avfd_cryo1["xmin"] < sr->slc[0].vertex.x) && (sr->slc[0].vertex.x < avfd_cryo1["xmax"]);
  		bool y = (avfd_cryo1["ymin"] < sr->slc[0].vertex.y) && (sr->slc[0].vertex.y < avfd_cryo1["ymax"]);
  		bool z = (avfd_cryo1["zmin"] < sr->slc[0].vertex.z) && (sr->slc[0].vertex.z < avfd_cryo1["zmax"]);
  		return(x && y && z);
  	}
  	);

  const Cut kActiveVolumeFDCryo2(
  	[](const caf::SRProxy* sr){
  		bool x = (avfd_cryo2["xmin"] < sr->slc[0].vertex.x) && (sr->slc[0].vertex.x < avfd_cryo2["xmax"]);
  		bool y = (avfd_cryo2["ymin"] < sr->slc[0].vertex.y) && (sr->slc[0].vertex.y < avfd_cryo2["ymax"]);
  		bool z = (avfd_cryo2["zmin"] < sr->slc[0].vertex.z) && (sr->slc[0].vertex.z < avfd_cryo2["zmax"]);
  		return(x && y && z);
  	}
  	);


} // namespace
