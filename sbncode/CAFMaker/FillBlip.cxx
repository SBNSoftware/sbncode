#include "sbncode/CAFMaker/FillBlip.h"

namespace caf
{
    void FillBlip(   const std::vector<blip::Blip>& LarBlips, std::vector<caf::SRBlip>& CAF_Blips)
    {
        for(blip::Blip const& CurrentBlip: LarBlips)
        {
            caf::SRBlip NewBlip;
            NewBlip.ID = CurrentBlip.ID;
            NewBlip.isValid = CurrentBlip.isValid;
            NewBlip.cryostat = CurrentBlip.Cryostat;
            NewBlip.TPC = CurrentBlip.TPC;
            NewBlip.nPlanes = CurrentBlip.NPlanes;
            NewBlip.maxWireSpan = CurrentBlip.MaxWireSpan;
            NewBlip.timeTick = CurrentBlip.TimeTick;
            NewBlip.time = CurrentBlip.Time;
            NewBlip.charge = CurrentBlip.Charge;
            NewBlip.energy = CurrentBlip.Energy/1000.; //convert to GeV
            NewBlip.energyESTAR = CurrentBlip.EnergyESTAR/1000.; //convert to GeV
            NewBlip.energyPSTAR = CurrentBlip.EnergyPSTAR/1000.; //convert to GeV
            NewBlip.proxTrkDist = CurrentBlip.ProxTrkDist; 
            NewBlip.proxTrkID = CurrentBlip.ProxTrkID;
            NewBlip.inCylinder = CurrentBlip.inCylinder;
            NewBlip.position.SetXYZ(CurrentBlip.Position.X(), CurrentBlip.Position.Y(), CurrentBlip.Position.Z());
            NewBlip.sigmaYZ = CurrentBlip.SigmaYZ;
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
      CAF_Blip.truthBlip.cryostat =LarBlip.truth.Cryostat;
      CAF_Blip.truthBlip.TPC =LarBlip.truth.TPC;
      CAF_Blip.truthBlip.time =LarBlip.truth.Time;
      CAF_Blip.truthBlip.driftTime =LarBlip.truth.DriftTime;
      CAF_Blip.truthBlip.energy =LarBlip.truth.Energy;
      CAF_Blip.truthBlip.depElectrons =LarBlip.truth.DepElectrons;
      CAF_Blip.truthBlip.numElectrons =LarBlip.truth.NumElectrons;
      CAF_Blip.truthBlip.leadG4ID =LarBlip.truth.LeadG4ID;
      CAF_Blip.truthBlip.leadG4Index =LarBlip.truth.LeadG4Index;
      CAF_Blip.truthBlip.leadG4PDG =LarBlip.truth.LeadG4PDG;
      CAF_Blip.truthBlip.leadCharge =LarBlip.truth.LeadCharge;
      CAF_Blip.truthBlip.position.SetXYZ(LarBlip.truth.Position.X(), LarBlip.truth.Position.Y(), LarBlip.truth.Position.Z());
      CAF_Blip.truthBlip.energy = CAF_Blip.truthBlip.energy/1000.; //convert to GeV
    }

    void FillBlipRealtedHitCluster(blip::Blip& LarBlip, caf::SRBlip &CAF_Blip)
    {
        int NumPlanes = sizeof(LarBlip.clusters)/sizeof(LarBlip.clusters[0]);
        for(int iPlane=0; iPlane<NumPlanes; iPlane++)
        {
            CAF_Blip.clusters[iPlane].ID = LarBlip.clusters[iPlane].ID;
            CAF_Blip.clusters[iPlane].isValid = LarBlip.clusters[iPlane].isValid;
            CAF_Blip.clusters[iPlane].centerChan = LarBlip.clusters[iPlane].CenterChan;
            CAF_Blip.clusters[iPlane].centerWire = LarBlip.clusters[iPlane].CenterWire;
            CAF_Blip.clusters[iPlane].isTruthMatched = LarBlip.clusters[iPlane].isTruthMatched;
            CAF_Blip.clusters[iPlane].isMatched = LarBlip.clusters[iPlane].isMatched;
            CAF_Blip.clusters[iPlane].deadWireSep = LarBlip.clusters[iPlane].DeadWireSep;
            CAF_Blip.clusters[iPlane].cryostat = LarBlip.clusters[iPlane].Cryostat;
            CAF_Blip.clusters[iPlane].TPC = LarBlip.clusters[iPlane].TPC;
            CAF_Blip.clusters[iPlane].plane = LarBlip.clusters[iPlane].Plane;
            CAF_Blip.clusters[iPlane].nHits = LarBlip.clusters[iPlane].NHits;
            CAF_Blip.clusters[iPlane].nWires = LarBlip.clusters[iPlane].NWires;
            CAF_Blip.clusters[iPlane].ADCs = LarBlip.clusters[iPlane].ADCs;
            CAF_Blip.clusters[iPlane].amplitude = LarBlip.clusters[iPlane].Amplitude;
            CAF_Blip.clusters[iPlane].charge = LarBlip.clusters[iPlane].Charge;
            CAF_Blip.clusters[iPlane].sigmaCharge = LarBlip.clusters[iPlane].SigmaCharge;
            CAF_Blip.clusters[iPlane].timeTick = LarBlip.clusters[iPlane].TimeTick;
            CAF_Blip.clusters[iPlane].time = LarBlip.clusters[iPlane].Time;
            CAF_Blip.clusters[iPlane].startTime = LarBlip.clusters[iPlane].StartTime;
            CAF_Blip.clusters[iPlane].endTime = LarBlip.clusters[iPlane].EndTime;
            CAF_Blip.clusters[iPlane].timespan = LarBlip.clusters[iPlane].Timespan;
            CAF_Blip.clusters[iPlane].RMS = LarBlip.clusters[iPlane].RMS;
            CAF_Blip.clusters[iPlane].startWire = LarBlip.clusters[iPlane].StartWire;
            CAF_Blip.clusters[iPlane].endWire = LarBlip.clusters[iPlane].EndWire;
            CAF_Blip.clusters[iPlane].nPulseTrainHits = LarBlip.clusters[iPlane].NPulseTrainHits;
            CAF_Blip.clusters[iPlane].goodnessOfFit = LarBlip.clusters[iPlane].GoodnessOfFit;
            CAF_Blip.clusters[iPlane].blipID = LarBlip.clusters[iPlane].BlipID;
            CAF_Blip.clusters[iPlane].edepID = LarBlip.clusters[iPlane].EdepID;
            //These are sets that need to become vectors so we need to do some loops
            for(auto HitID : LarBlip.clusters[iPlane].HitIDs) CAF_Blip.clusters[iPlane].hitIDs.push_back(HitID);
            for(auto Wire : LarBlip.clusters[iPlane].Wires) CAF_Blip.clusters[iPlane].wires.push_back(Wire);
            for(auto Chan : LarBlip.clusters[iPlane].Chans) CAF_Blip.clusters[iPlane].chans.push_back(Chan);
            for(auto G4ID : LarBlip.clusters[iPlane].G4IDs) CAF_Blip.clusters[iPlane].G4IDs.push_back(G4ID);
        }
    }
}
