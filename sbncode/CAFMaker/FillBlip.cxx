#include "sbncode/CAFMaker/FillBlip.h"

namespace caf
{
    void FillBlip(   const std::vector<blip::Blip>& LarBlips, std::vector<caf::SRBlip>& CAF_Blips)
    {
        int NBlips = LarBlips.size();
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
            NewBlip.Position.SetXYZ(CurrentBlip.Position.X(), CurrentBlip.Position.Y(), CurrentBlip.Position.Z());
            NewBlip.SigmaYZ = CurrentBlip.SigmaYZ;
            NewBlip.dX = CurrentBlip.dX;
            NewBlip.dYZ = CurrentBlip.dYZ;
            if(CurrentBlip.truth.ID >= 0 ) //MC Blip
            {
                FillMCTruthBlip( CurrentBlip, NewBlip );
            }
            FillBlipRealtedHitCluster( CurrentBlip, NewBlip );
            CAF_Blips.push_back(NewBlip);
        }

    }

    void FillMCTruthBlip( blip::Blip& LarBlip, caf::SRBlip &CAF_Blip )
    {
      CAF_Blip.truthBlip.ID = LarBlip.truth.ID;
      CAF_Blip.truthBlip.Cryostat =LarBlip.truth.Cryostat;
      CAF_Blip.truthBlip.TPC =LarBlip.truth.TPC;
      CAF_Blip.truthBlip.Time =LarBlip.truth.Time;
      CAF_Blip.truthBlip.DriftTime =LarBlip.truth.DriftTime;
      CAF_Blip.truthBlip.Energy =LarBlip.truth.Energy;
      CAF_Blip.truthBlip.DepElectrons =LarBlip.truth.DepElectrons;
      CAF_Blip.truthBlip.NumElectrons =LarBlip.truth.NumElectrons;
      CAF_Blip.truthBlip.LeadG4ID =LarBlip.truth.LeadG4ID;
      CAF_Blip.truthBlip.LeadG4Index =LarBlip.truth.LeadG4Index;
      CAF_Blip.truthBlip.LeadG4PDG =LarBlip.truth.LeadG4PDG;
      CAF_Blip.truthBlip.LeadCharge =LarBlip.truth.LeadCharge;
      CAF_Blip.truthBlip.Position.SetXYZ(LarBlip.truth.Position.X(), LarBlip.truth.Position.Y(), LarBlip.truth.Position.Z());
      CAF_Blip.truthBlip.Energy = CAF_Blip.truthBlip.Energy/1000.; //convert to GeV
    }

    void FillBlipRealtedHitCluster(blip::Blip& LarBlip, caf::SRBlip &CAF_Blip)
    {
        int NumPlanes = sizeof(LarBlip.clusters)/sizeof(LarBlip.clusters[0]);
        for(int iPlane=0; iPlane<NumPlanes; iPlane++)
        {
            CAF_Blip.clusters[iPlane].ID = LarBlip.clusters[iPlane].ID;
            CAF_Blip.clusters[iPlane].isValid = LarBlip.clusters[iPlane].isValid;
            CAF_Blip.clusters[iPlane].CenterChan = LarBlip.clusters[iPlane].CenterChan;
            CAF_Blip.clusters[iPlane].CenterWire = LarBlip.clusters[iPlane].CenterWire;
            CAF_Blip.clusters[iPlane].isTruthMatched = LarBlip.clusters[iPlane].isTruthMatched;
            CAF_Blip.clusters[iPlane].isMatched = LarBlip.clusters[iPlane].isMatched;
            CAF_Blip.clusters[iPlane].DeadWireSep = LarBlip.clusters[iPlane].DeadWireSep;
            CAF_Blip.clusters[iPlane].Cryostat = LarBlip.clusters[iPlane].Cryostat;
            CAF_Blip.clusters[iPlane].TPC = LarBlip.clusters[iPlane].TPC;
            CAF_Blip.clusters[iPlane].Plane = LarBlip.clusters[iPlane].Plane;
            CAF_Blip.clusters[iPlane].NHits = LarBlip.clusters[iPlane].NHits;
            CAF_Blip.clusters[iPlane].NWires = LarBlip.clusters[iPlane].NWires;
            CAF_Blip.clusters[iPlane].ADCs = LarBlip.clusters[iPlane].ADCs;
            CAF_Blip.clusters[iPlane].Amplitude = LarBlip.clusters[iPlane].Amplitude;
            CAF_Blip.clusters[iPlane].Charge = LarBlip.clusters[iPlane].Charge;
            CAF_Blip.clusters[iPlane].SigmaCharge = LarBlip.clusters[iPlane].SigmaCharge;
            CAF_Blip.clusters[iPlane].TimeTick = LarBlip.clusters[iPlane].TimeTick;
            CAF_Blip.clusters[iPlane].Time = LarBlip.clusters[iPlane].Time;
            CAF_Blip.clusters[iPlane].StartTime = LarBlip.clusters[iPlane].StartTime;
            CAF_Blip.clusters[iPlane].EndTime = LarBlip.clusters[iPlane].EndTime;
            CAF_Blip.clusters[iPlane].Timespan = LarBlip.clusters[iPlane].Timespan;
            CAF_Blip.clusters[iPlane].RMS = LarBlip.clusters[iPlane].RMS;
            CAF_Blip.clusters[iPlane].StartWire = LarBlip.clusters[iPlane].StartWire;
            CAF_Blip.clusters[iPlane].EndWire = LarBlip.clusters[iPlane].EndWire;
            CAF_Blip.clusters[iPlane].NPulseTrainHits = LarBlip.clusters[iPlane].NPulseTrainHits;
            CAF_Blip.clusters[iPlane].GoodnessOfFit = LarBlip.clusters[iPlane].GoodnessOfFit;
            CAF_Blip.clusters[iPlane].BlipID = LarBlip.clusters[iPlane].BlipID;
            CAF_Blip.clusters[iPlane].EdepID = LarBlip.clusters[iPlane].EdepID;
            //These are sets that need to become vectors so we need to do some loops
            for(auto HitID : LarBlip.clusters[iPlane].HitIDs) CAF_Blip.clusters[iPlane].HitIDs.push_back(HitID);
            for(auto Wire : LarBlip.clusters[iPlane].Wires) CAF_Blip.clusters[iPlane].Wires.push_back(Wire);
            for(auto Chan : LarBlip.clusters[iPlane].Chans) CAF_Blip.clusters[iPlane].Chans.push_back(Chan);
            for(auto G4ID : LarBlip.clusters[iPlane].G4IDs) CAF_Blip.clusters[iPlane].G4IDs.push_back(G4ID);
        }
    }
}
