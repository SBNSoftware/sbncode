#include "sbncode/CAFMaker/FillBlip.h"

namespace caf
{
    void FillBlip(   const std::vector<blip::Blip>& LarBlips, std::vector<caf::SRBlip>& CAF_Blips)
    {
        int NBlips = LarBlips->size();
        for(int iBlip=0; iBlip<NBlips; iBlip++)
        {
            blip::Blip CurrentBlip = LarBlips[iBlip];
            caf::SRBlip NewBlip;
            NewBlip.ID = CurrentBlip.ID;
            NewBlip.isValid = CurrentBlip.isValid;
            NewBlip.Cryostat = CurrentBlip.Cryostat;
            NewBlip.TPC = CurrentBlip.TPC;
            NewBlip.NPlanes = CurrentBlip.NPlanes;
            NewBlip.MaxWireSpan = CurrentBlip.MaxWireSpan;
            NewBlip.TimeTick = CurrentBlip.TimeTick;
            NewBlip.Time = CurrentBlip.Time;
            NewBlip.Charge = CurrentBlip.Charge;
            NewBlip.Energy = CurrentBlip.Energy/1000.; //convert to GeV
            NewBlip.EnergyESTAR = CurrentBlip.EnergyESTAR/1000.; //convert to GeV
            NewBlip.EnergyPSTAR = CurrentBlip.EnergyPSTAR/1000.; //convert to GeV
            NewBlip.ProxTrkDist = CurrentBlip.ProxTrkDist; 
            NewBlip.ProxTrkID = CurrentBlip.ProxTrkID;
            NewBlip.inCylinder = CurrentBlip.inCylinder;
            NewBlip.Position = CurrentBlip.Position;
            NewBlip.SigmaYZ = CurrentBlip.SigmaYZ;
            NewBlip.dX = CurrentBlip.dX;
            NewBlip.dYZ = CurrentBlip.dYZ;
            if(CurrentBlip.truthBlip.ID >= 0 ) //MC Blip
            {
                FillMCTruthBlip( &CurrentBlip, &NewBlip );
            }
            FillBlipRealtedHitCluster( &CurrentBlip, &NewBlip );
        }

    }

    void FillMCTruthBlip( blip::Blip& LarBlip, caf::SRBlip &CAF_Blip )
    {
        CAF_Blip->truthBlip = LarBlip->truthBlip;
        CAF_Blip->truthBlip.Energy = CAF_Blip->truthBlip.Energy/1000.; //convert to GeV
    }

    void FillBlipRealtedHitCluster(blip::Blip& LarBlip, caf::SRBlip &CAF_Blip)
    {
        int NumPlanes = sizeof(LarBlip->clusters)/sizeof(LarBlip->clusters[0]);
        for(int iPlane=0; iPlane<NumPlanes; iPlane++)
        {
            CAF_Blip->clusters[iPlane] = LarBlip->clusters[iPlane]
        }
    }
}