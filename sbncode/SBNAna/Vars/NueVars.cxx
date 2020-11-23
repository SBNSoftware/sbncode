#include "SBNAna/Vars/NueVars.h"
#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{
  // Get the Largest shower in the slice
  const Var kLargestRecoShowerIdx(
      [](const caf::SRSliceProxy *slc)
      {
        int bestIdx(-1);
        double maxEnergy(-1);

        for (unsigned int i=0; i< slc->reco.nshw; i++) {
          auto const& shw = slc->reco.shw[i];
          if (shw.bestplane_energy > maxEnergy){
            bestIdx = i;
            maxEnergy = shw.bestplane_energy;
          }
        }
        return bestIdx;
      });

  // Currently assumes shw 0 is the primary
  const Var kRecoShower_BestEnergy(
			[](const caf::SRSliceProxy *slc)
			{
			  double energy = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    energy = slc->reco.shw[largestShwIdx].bestplane_energy;
			  }
			  return energy;
			});

  // Currently assumes shw 0 is the primary
  const Var kRecoShower_BestdEdx(
			[](const caf::SRSliceProxy *slc)
			{
			  double dedx = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    dedx = slc->reco.shw[largestShwIdx].bestplane_dEdx;
			  }
			  return dedx;
			});

  const Var kRecoShower_ConversionGap(
			[](const caf::SRSliceProxy *slc)
			{
        double gap = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
        if ( largestShwIdx!=-1 ){
          gap = slc->reco.shw[largestShwIdx].conversion_gap;
        }
			  return gap;
			});

  const Var kRecoShower_Density(
			[](const caf::SRSliceProxy *slc)
			{
			  double density = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
        if ( largestShwIdx!=-1 ){
			    density = slc->reco.shw[largestShwIdx].density;
			  }
			  return density;
			});

  const Var kRecoShower_Energy(
			[](const caf::SRSliceProxy *slc)
			{
			  double energy = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    energy = slc->reco.shw[largestShwIdx].energy_plane1; // so far taking whatever plane 1 is and first shw
			    // energy = slc->reco.shw[largestShwIdx].energy.x; // so far taking whatever plane 1 is and first shw
			  }
			  return energy;
			});

  const Var kRecoShower_Length(
			       [](const caf::SRSliceProxy *slc)
			       {
				 double len = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
				   len = slc->reco.shw[largestShwIdx].len;
				 }
				 return len;
			       });

    const Var kRecoShower_OpenAngle(
			[](const caf::SRSliceProxy *slc)
			{
			  double openangle = -5.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    openangle = slc->reco.shw[largestShwIdx].open_angle;
			  }
			  return openangle;
			});


    const Var kRecoShower_StartX(
			[](const caf::SRSliceProxy *slc)
			{
			  double vtx = -9999.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    vtx = slc->reco.shw[largestShwIdx].start.x;
			  }
			  return vtx;
			});

    const Var kRecoShower_StartY(
			[](const caf::SRSliceProxy *slc)
			{
			  double vtx = -9999.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    vtx = slc->reco.shw[largestShwIdx].start.y;
			  }
			  return vtx;
			});

    const Var kRecoShower_StartZ(
			[](const caf::SRSliceProxy *slc)
			{
			  double vtx = -9999.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    vtx = slc->reco.shw[largestShwIdx].start.z;
			  }
			  return vtx;
			});

    const Var kRecoShower_EndX(
			[](const caf::SRSliceProxy *slc)
			{
			  double end = -9999.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    double start = slc->reco.shw[largestShwIdx].start.x;
			    double dir = slc->reco.shw[largestShwIdx].dir.x;
			    double len = slc->reco.shw[largestShwIdx].len;
			    end = start + (dir*len);
			  }
			  return end;
			});

    const Var kRecoShower_EndY(
			[](const caf::SRSliceProxy *slc)
			{
			  double end = -9999.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    double start = slc->reco.shw[largestShwIdx].start.y;
			    double dir = slc->reco.shw[largestShwIdx].dir.y;
			    double len = slc->reco.shw[largestShwIdx].len;
			    end = start + (dir*len);
			  }
			  return end;
			});

    const Var kRecoShower_EndZ(
			[](const caf::SRSliceProxy *slc)
			{
			  double end = -9999.0;
        const int largestShwIdx(kLargestRecoShowerIdx(slc));
			  if ( largestShwIdx!=-1 ){
			    double start = slc->reco.shw[largestShwIdx].start.z;
			    double dir = slc->reco.shw[largestShwIdx].dir.z;
			    double len = slc->reco.shw[largestShwIdx].len;
			    end = start + (dir*len);
			  }
			  return end;
			});

}
