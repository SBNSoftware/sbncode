#include "SBNAna/Vars/NueVars.h"
#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{

  // Currently assumes shw 0 is the primary
  const Var kRecoShower_BestEnergy(
			[](const caf::SRSliceProxy *slc)
			{
			  double energy = -5.0;
			  if ( slc->reco.nshw ){
			    energy = slc->reco.shw[0].bestplane_energy;
			  }
			  return energy;
			});

  const Var kRecoShower_ConversionGap(
			[](const caf::SRSliceProxy *slc)
			{
                          double gap = -5.0;
                          if ( slc->reco.nshw ){
                            double x = slc->vertex.x - slc->reco.shw[0].start.x;
                            double y = slc->vertex.y - slc->reco.shw[0].start.y;
                            double z = slc->vertex.z - slc->reco.shw[0].start.z;
                            gap  = sqrt(x*x + y*y + z*z);	
                          }

			  return gap;
			});

  const Var kRecoShower_Density(
			[](const caf::SRSliceProxy *slc)
			{
			  double density = -5.0;
			  if ( slc->reco.nshw && slc->reco.shw[0].len >0 ){
			    density = slc->reco.shw[0].bestplane_energy / slc->reco.shw[0].len;
			  }
			  return density;
			});

  const Var kRecoShower_Energy(
			[](const caf::SRSliceProxy *slc)
			{
			  double energy = -5.0;
			  if ( slc->reco.nshw ){
			    energy = slc->reco.shw[0].energy[0]; // so far taking whatever plane 1 is
			  }
			  return energy;
			});

  const Var kRecoShower_Length(
			       [](const caf::SRSliceProxy *slc)
			       {
				 double len = -5.0;
				 if ( slc->reco.nshw ){
				   len = slc->reco.shw[0].len;
				 }
				 return len;
			       });

    const Var kRecoShower_OpenAngle(
			[](const caf::SRSliceProxy *slc)
			{
			  double openangle = -5.0;
			  if ( slc->reco.nshw ){
			    openangle = slc->reco.shw[0].open_angle;
			  }
			  return openangle;
			});


    const Var kRecoShower_StartX(
			[](const caf::SRSliceProxy *slc)
			{
			  double vtx = -5.0;
			  if ( slc->reco.nshw ){
			    vtx = slc->reco.shw[0].start.x;
			  }
			  return vtx;
			});

    const Var kRecoShower_StartY(
			[](const caf::SRSliceProxy *slc)
			{
			  double vtx = -5.0;
			  if ( slc->reco.nshw ){
			    vtx = slc->reco.shw[0].start.y;
			  }
			  return vtx;
			});

    const Var kRecoShower_StartZ(
			[](const caf::SRSliceProxy *slc)
			{
			  double vtx = -5.0;
			  if ( slc->reco.nshw ){
			    vtx = slc->reco.shw[0].start.z;
			  }
			  return vtx;
			});

}
