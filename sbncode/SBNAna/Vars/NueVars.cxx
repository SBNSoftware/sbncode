#include "SBNAna/Vars/NueVars.h"
#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{

  // Currently assumes shw 0 is the primary
  const Var kRecoShower_BestEnergy(
			[](const caf::SRProxy *sr)
			{
			  double energy = -5.0;
			  if ( sr->reco.nshw ){
			    energy = sr->reco.shw[0].bestplane_energy;
			  }
			  return energy;
			});

  // Placeholder. Need better definition which needs better caf filling.
  const Var kRecoShower_ConversionGap(
			[](const caf::SRProxy *sr)
			{
				double gap = -5.0;
				if ( sr->reco.nshw ){
					double x = sr->slc.vertex.x - sr->reco.shw[0].start.x;
					double y = sr->slc.vertex.y - sr->reco.shw[0].start.y;
					double z = sr->slc.vertex.z - sr->reco.shw[0].start.z;
					gap  = sqrt(x*x + y*y + z*z);	
				}
			  return gap;
			});

  const Var kRecoShower_Density(
			[](const caf::SRProxy *sr)
			{
			  double density = -5.0;
			  if ( sr->reco.nshw && sr->reco.shw[0].len >0 ){
			    density = sr->reco.shw[0].bestplane_energy / sr->reco.shw[0].len;
			  }
			  return density;
			});

  const Var kRecoShower_Energy(
			[](const caf::SRProxy *sr)
			{
			  double energy = -5.0;
			  if ( sr->reco.nshw ){
			    energy = sr->reco.shw[0].energy[0]; // so far taking whatever plane 1 is
			  }
			  return energy;
			});

  const Var kRecoShower_Length(
			       [](const caf::SRProxy *sr)
			       {
				 double len = -5.0;
				 if ( sr->reco.nshw ){
				   len = sr->reco.shw[0].len;
				 }
				 return len;
			       });

    const Var kRecoShower_OpenAngle(
			[](const caf::SRProxy *sr)
			{
			  double openangle = -5.0;
			  if ( sr->reco.nshw ){
			    openangle = sr->reco.shw[0].open_angle;
			  }
			  return openangle;
			});


    const Var kRecoShower_StartX(
			[](const caf::SRProxy *sr)
			{
			  double vtx = -5.0;
			  if ( sr->reco.nshw ){
			    vtx = sr->reco.shw[0].start.x;
			  }
			  return vtx;
			});

    const Var kRecoShower_StartY(
			[](const caf::SRProxy *sr)
			{
			  double vtx = -5.0;
			  if ( sr->reco.nshw ){
			    vtx = sr->reco.shw[0].start.y;
			  }
			  return vtx;
			});

    const Var kRecoShower_StartZ(
			[](const caf::SRProxy *sr)
			{
			  double vtx = -5.0;
			  if ( sr->reco.nshw ){
			    vtx = sr->reco.shw[0].start.z;
			  }
			  return vtx;
			});

}
