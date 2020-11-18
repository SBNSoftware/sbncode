#include "SBNAna/Cuts/NueCuts.h"
#include "SBNAna/Vars/NueVars.h"

#include "StandardRecord/Proxy/SRProxy.h"

// icarus active volume cryo 1
std::map<std::string,double> avfd =
{
	{"xmin", -368.49},
	{"xmax", -71.94},
	{"ymin", -181.86},
	{"ymax", 134.96},
	{"zmin", -894.95},
	{"zmax", 894.95}
};

namespace ana{

	// Basic reconstruction 
	const Cut kRecoShower(
		[](const caf::SRSliceProxy* slc)
		{
		  return (slc->reco.nshw > 0 && // need a shower
				  slc->reco.shw[0].bestplane_energy > 0 && // nothing is terribly wrong
				  slc->reco.shw[0].len > 0 );
	}
	);

	// Basic reconstruction 
	const Cut kNueBasicCut(
		[](const caf::SRSliceProxy* slc)
		{
		  return (slc->reco.shw[0].bestplane_energy < 250. &&
				  slc->reco.shw[0].bestplane_dEdx < 2.7 &&
				  slc->reco.shw[0].len < 42.);
	}
	);

	// // Cut currently not working as it wants a caf::SRProxy and not caf::SRSliceProxy
	// // Workaround is definiting it with the explicit branch i.e. slc->reco.shw[0].start.z instead of kRecoShower_StartZ
	// const Cut kNueContainedFD(
	// 	[](const caf::SRSliceProxy* slc){

	// 		bool xstart = (avfd_cryo1["xmin"] < kRecoShower_StartX) && (kRecoShower_StartX < avfd_cryo1["xmax"]);
	// 		bool xend   = (avfd_cryo1["xmin"] < kRecoShower_EndX) && (kRecoShower_EndX < avfd_cryo1["xmax"]);

	// 		bool ystart = (avfd_cryo1["ymin"] < kRecoShower_StartY) && (kRecoShower_StartY < avfd_cryo1["ymax"]);
	// 		bool yend   = (avfd_cryo1["ymin"] < kRecoShower_EndY) && (kRecoShower_EndY < avfd_cryo1["ymax"]);

	// 		bool zstart = (avfd_cryo1["zmin"] < kRecoShower_StartZ) && (kRecoShower_StartZ < avfd_cryo1["zmax"]);
	// 		bool zend   = (avfd_cryo1["zmin"] < kRecoShower_EndZ) && (kRecoShower_EndZ < avfd_cryo1["zmax"]);

	// 		return (xstart && xend && ystart && yend && zstart && zend);
 //    }
 //    );

	const Cut kNueContainedFD(
		[](const caf::SRSliceProxy* slc){

			double this_endx = -9999.0;
			double this_endy = -9999.0;
			double this_endz = -9999.0;
			if ( slc->reco.nshw ){
			    this_endx = slc->reco.shw[0].start.x + (slc->reco.shw[0].dir.x * slc->reco.shw[0].len);
			    this_endy = slc->reco.shw[0].start.y + (slc->reco.shw[0].dir.y * slc->reco.shw[0].len);
			    this_endz = slc->reco.shw[0].start.z + (slc->reco.shw[0].dir.z * slc->reco.shw[0].len);
			}

			bool startx = (avfd["xmin"] < slc->reco.shw[0].start.x) && (slc->reco.shw[0].start.x < avfd["xmax"]);
			bool endx   = (avfd["xmin"] < this_endx) && (this_endx < avfd["xmax"]);

			bool starty = (avfd["ymin"] < slc->reco.shw[0].start.y) && (slc->reco.shw[0].start.y < avfd["ymax"]);
			bool endy   = (avfd["ymin"] < this_endy) && (this_endy < avfd["ymax"]);

			bool startz = (avfd["zmin"] < slc->reco.shw[0].start.z) && (slc->reco.shw[0].start.z < avfd["zmax"]);
			bool endz   = (avfd["zmin"] < this_endz) && (this_endz < avfd["zmax"]);

			return (startx && endx && starty && endy && startz && endz);
    }
    );

}
