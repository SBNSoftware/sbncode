#include "SBNAna/Cuts/NueCuts.h"
#include "SBNAna/Vars/NueVars.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include "SBNAna/Cuts/VolumeDefinitions.h"

namespace ana{

  // Basic reconstruction
  const Cut kRecoShower(
    [](const caf::SRSliceProxy* slc)
    {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      return ( slc->reco.shw[largestShwIdx].bestplane_energy > 0 && // nothing is terribly wrong
          slc->reco.shw[largestShwIdx].bestplane_dEdx > 0 &&
          slc->reco.shw[largestShwIdx].conversion_gap > 0 );
  }
  );

  // Basic reconstruction
  const Cut kNueBasicCut(
    [](const caf::SRSliceProxy* slc)
    {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      return (slc->reco.shw[largestShwIdx].bestplane_energy > 200. &&
          slc->reco.shw[largestShwIdx].bestplane_dEdx < 3 &&
          slc->reco.shw[largestShwIdx].conversion_gap < 3 );
  }
  );

  const Cut kShowerEnergyCut(
    [](const caf::SRSliceProxy* slc)
    {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      return (slc->reco.shw[largestShwIdx].bestplane_energy > 200.f );
    }
    );

  const Cut kShowerdEdxCut(
    [](const caf::SRSliceProxy* slc)
    {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      return (slc->reco.shw[largestShwIdx].bestplane_dEdx < 3.f);
    }
    );

  const Cut kShowerConvGapCut(
    [](const caf::SRSliceProxy* slc)
    {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      return (slc->reco.shw[largestShwIdx].conversion_gap > 0.f && slc->reco.shw[largestShwIdx].conversion_gap < 3.f);
    }
    );

  // Basic reconstruction
  const Cut kNueNumShowersCut(
    [](const caf::SRSliceProxy* slc)
    {
      return ((unsigned int)kRecoShowers_EnergyCut(slc) > 0);
    }
    );

  const Cut kNueTrackLenCut(
    [](const caf::SRSliceProxy* slc)
    {
      return kLongestTrackLength(slc) < 50.;
    }
    );

  const Cut kShowerDensityCut(
    [](const caf::SRSliceProxy* slc)
    {
      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      return (slc->reco.shw[largestShwIdx].density > 6);
    }
    );

  // // Cut currently not working as it wants a caf::SRProxy and not caf::SRSliceProxy
  // // Workaround is definiting it with the explicit branch i.e. slc->reco.shw[0].start.z instead of kRecoShower_StartZ
  // const Cut kNueContainedFD(
  //  [](const caf::SRSliceProxy* slc){

  //    bool xstart = (avfd_cryo1_cryo1.xmin < kRecoShower_StartX) && (kRecoShower_StartX < avfd_cryo1_cryo1.xmax);
  //    bool xend   = (avfd_cryo1_cryo1.xmin < kRecoShower_EndX) && (kRecoShower_EndX < avfd_cryo1_cryo1.xmax);

  //    bool ystart = (avfd_cryo1_cryo1.ymin < kRecoShower_StartY) && (kRecoShower_StartY < avfd_cryo1_cryo1.ymax);
  //    bool yend   = (avfd_cryo1_cryo1.ymin < kRecoShower_EndY) && (kRecoShower_EndY < avfd_cryo1_cryo1.ymax);

  //    bool zstart = (avfd_cryo1_cryo1.zmin < kRecoShower_StartZ) && (kRecoShower_StartZ < avfd_cryo1_cryo1.zmax);
  //    bool zend   = (avfd_cryo1_cryo1.zmin < kRecoShower_EndZ) && (kRecoShower_EndZ < avfd_cryo1_cryo1.zmax);

  //    return (xstart && xend && ystart && yend && zstart && zend);
 //    }
 //    );

  // TODO: find a better way to set the AV so we do not need to replicate code
  // shw.end ha been added to the CAF so should be trivial for future iterations
  const Cut kNueContainedND(
    [](const caf::SRSliceProxy* slc){

      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      double this_endx = slc->reco.shw[largestShwIdx].start.x + (slc->reco.shw[largestShwIdx].dir.x * slc->reco.shw[largestShwIdx].len);
      double this_endy = slc->reco.shw[largestShwIdx].start.y + (slc->reco.shw[largestShwIdx].dir.y * slc->reco.shw[largestShwIdx].len);
      double this_endz = slc->reco.shw[largestShwIdx].start.z + (slc->reco.shw[largestShwIdx].dir.z * slc->reco.shw[largestShwIdx].len);

      bool startx = (fvndAbs.xmin < std::abs(slc->reco.shw[largestShwIdx].start.x)) && (std::abs(slc->reco.shw[largestShwIdx].start.x) < fvndAbs.xmax);
      bool endx   = (fvndAbs.xmin < std::abs(this_endx)) && (std::abs(this_endx) < fvndAbs.xmax);

      bool starty = (fvndAbs.ymin < slc->reco.shw[largestShwIdx].start.y) && (slc->reco.shw[largestShwIdx].start.y < fvndAbs.ymax);
      bool endy   = (fvndAbs.ymin < this_endy) && (this_endy < fvndAbs.ymax);

      bool startz = (fvndAbs.zmin < slc->reco.shw[largestShwIdx].start.z) && (slc->reco.shw[largestShwIdx].start.z < fvndAbs.zmax);
      bool endz   = (fvndAbs.zmin < this_endz) && (this_endz < fvndAbs.zmax);

      return (startx && endx && starty && endy && startz && endz);
    }
    );

  const Cut kNueContainedFD(
    [](const caf::SRSliceProxy* slc){

      const int largestShwIdx(kLargestRecoShowerIdx(slc));
      if ( largestShwIdx==-1 )
        return false;

      double this_endx = slc->reco.shw[largestShwIdx].start.x + (slc->reco.shw[largestShwIdx].dir.x * slc->reco.shw[largestShwIdx].len);
      double this_endy = slc->reco.shw[largestShwIdx].start.y + (slc->reco.shw[largestShwIdx].dir.y * slc->reco.shw[largestShwIdx].len);
      double this_endz = slc->reco.shw[largestShwIdx].start.z + (slc->reco.shw[largestShwIdx].dir.z * slc->reco.shw[largestShwIdx].len);


      bool startx = (fvfd_cryo1.xmin < slc->reco.shw[largestShwIdx].start.x) && (slc->reco.shw[largestShwIdx].start.x < fvfd_cryo1.xmax);
      bool endx   = (fvfd_cryo1.xmin < this_endx) && (this_endx < fvfd_cryo1.xmax);

      bool starty = (fvfd_cryo1.ymin < slc->reco.shw[largestShwIdx].start.y) && (slc->reco.shw[largestShwIdx].start.y < fvfd_cryo1.ymax);
      bool endy   = (fvfd_cryo1.ymin < this_endy) && (this_endy < fvfd_cryo1.ymax);

      bool startz = (fvfd_cryo1.zmin < slc->reco.shw[largestShwIdx].start.z) && (slc->reco.shw[largestShwIdx].start.z < fvfd_cryo1.zmax);
      bool endz   = (fvfd_cryo1.zmin < this_endz) && (this_endz < fvfd_cryo1.zmax);

      return (startx && endx && starty && endy && startz && endz);
    }
    );

}
