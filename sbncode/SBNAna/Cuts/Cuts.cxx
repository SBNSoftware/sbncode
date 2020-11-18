#include "SBNAna/Cuts/Cuts.h"
#include "SBNAna/Vars/Vars.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include <map>
#include <string>


namespace ana{

  // Dummy cut example for testing
  const SpillCut kFirstEvents = kEvt < 10;

  const SpillCut kFlashTrigger(
    [](const caf::SRSpillProxy* sr){
      return ( sr->pass_flashtrig == 1);
    }
    );

  const Cut kActiveVolumeND(
  	[](const caf::SRSliceProxy* slc){
  		bool x = (avnd["xmin"] < abs(slc->vertex.x)) && (abs(slc->vertex.x) < avnd["xmax"]);
  		bool y = (avnd["ymin"] < slc->vertex.y) && (slc->vertex.y < avnd["ymax"]);
  		bool z = (avnd["zmin"] < slc->vertex.z) && (slc->vertex.z < avnd["zmax"]);
  		return(x && y && z);
  	}
  	);

  const Cut kActiveVolumeFDCryo1(
  	[](const caf::SRSliceProxy* slc){
  		bool x = (avfd_cryo1["xmin"] < slc->vertex.x) && (slc->vertex.x < avfd_cryo1["xmax"]);
  		bool y = (avfd_cryo1["ymin"] < slc->vertex.y) && (slc->vertex.y < avfd_cryo1["ymax"]);
  		bool z = (avfd_cryo1["zmin"] < slc->vertex.z) && (slc->vertex.z < avfd_cryo1["zmax"]);
  		return(x && y && z);
  	}
  	);

  const Cut kActiveVolumeFDCryo2(
        [](const caf::SRSliceProxy* slc){
  		bool x = (avfd_cryo2["xmin"] < slc->vertex.x) && (slc->vertex.x < avfd_cryo2["xmax"]);
  		bool y = (avfd_cryo2["ymin"] < slc->vertex.y) && (slc->vertex.y < avfd_cryo2["ymax"]);
  		bool z = (avfd_cryo2["zmin"] < slc->vertex.z) && (slc->vertex.z < avfd_cryo2["zmax"]);
  		return(x && y && z);
  	}
  	);

} // namespace
