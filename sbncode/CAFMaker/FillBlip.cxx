#include "sbncode/CAFMaker/FillBlip.h"

namespace caf
{
    void FillBlip(   const std::vector<blip::Blip>& LAr_Blips, std::vector<caf::SRBlip>& CAF_Blips)
    {
        for(blip::Blip const& CurrentBlip: LAr_Blips)
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
                FillMCTruthBlip( CurrentBlip.truth, NewBlip.truthBlip );
            }
            for(int iPlane=0; iPlane<int(CAF_Blip.kNplanes); iPlane++)
            {
                FillBlipRealtedHitCluster( CurrentBlip.clusters[iPlane], NewBlip.clusters[iPlane] );
            }
            CAF_Blips.push_back(NewBlip);
        }

    }

    void FillMCTruthBlip( blip::TrueBlip const & TrueLAr_Blip, caf::SRBlipTrueBlip &TrueCAF_Blip )
    {
      TrueCAF_Blip.ID = TrueLAr_Blip.ID;
      TrueCAF_Blip.cryostat = TrueLAr_Blip.Cryostat;
      TrueCAF_Blip.TPC = TrueLAr_Blip.TPC;
      TrueCAF_Blip.time = TrueLAr_Blip.Time;
      TrueCAF_Blip.driftTime = TrueLAr_Blip.DriftTime;
      TrueCAF_Blip.energy = TrueLAr_Blip.Energy;
      TrueCAF_Blip.depElectrons = TrueLAr_Blip.DepElectrons;
      TrueCAF_Blip.numElectrons = TrueLAr_Blip.NumElectrons;
      TrueCAF_Blip.leadG4ID = TrueLAr_Blip.LeadG4ID;
      TrueCAF_Blip.leadG4Index = TrueLAr_Blip.LeadG4Index;
      TrueCAF_Blip.leadG4PDG = TrueLAr_Blip.LeadG4PDG;
      TrueCAF_Blip.leadCharge = TrueLAr_Blip.LeadCharge;
      TrueCAF_Blip.position.SetXYZ(TrueLAr_Blip.Position.X(), TrueLAr_Blip.Position.Y(), TrueLAr_Blip.Position.Z());
      TrueCAF_Blip.energy = TrueCAF_Blip.energy/1000.; //convert to GeV
    }

    void FillBlipRealtedHitCluster(blip::HitClust const & LAr_HitClust, caf::SRBlipHitClust &CAF_HitClust)
    {
            CAF_HitClust.ID = LAr_HitClust.ID;
            CAF_HitClust.isValid = LAr_HitClust.isValid;
            CAF_HitClust.centerChan = LAr_HitClust.CenterChan;
            CAF_HitClust.centerWire = LAr_HitClust.CenterWire;
            CAF_HitClust.isTruthMatched = LAr_HitClust.isTruthMatched;
            CAF_HitClust.isMatched = LAr_HitClust.isMatched;
            CAF_HitClust.deadWireSep = LAr_HitClust.DeadWireSep;
            CAF_HitClust.cryostat = LAr_HitClust.Cryostat;
            CAF_HitClust.TPC = LAr_HitClust.TPC;
            CAF_HitClust.plane = LAr_HitClust.Plane;
            CAF_HitClust.nHits = LAr_HitClust.NHits;
            CAF_HitClust.nWires = LAr_HitClust.NWires;
            CAF_HitClust.ADCs = LAr_HitClust.ADCs;
            CAF_HitClust.amplitude = LAr_HitClust.Amplitude;
            CAF_HitClust.charge = LAr_HitClust.Charge;
            CAF_HitClust.sigmaCharge = LAr_HitClust.SigmaCharge;
            CAF_HitClust.timeTick = LAr_HitClust.TimeTick;
            CAF_HitClust.time = LAr_HitClust.Time;
            CAF_HitClust.startTime = LAr_HitClust.StartTime;
            CAF_HitClust.endTime = LAr_HitClust.EndTime;
            CAF_HitClust.timespan = LAr_HitClust.Timespan;
            CAF_HitClust.RMS = LAr_HitClust.RMS;
            CAF_HitClust.startWire = LAr_HitClust.StartWire;
            CAF_HitClust.endWire = LAr_HitClust.EndWire;
            CAF_HitClust.nPulseTrainHits = LAr_HitClust.NPulseTrainHits;
            CAF_HitClust.goodnessOfFit = LAr_HitClust.GoodnessOfFit;
            CAF_HitClust.blipID = LAr_HitClust.BlipID;
            CAF_HitClust.edepID = LAr_HitClust.EdepID;
            //These are sets that need to become vectors so we need to do some loops
            for(auto HitID : LAr_HitClust.HitIDs) CAF_HitClust.hitIDs.push_back(HitID);
            for(auto Wire : LAr_HitClust.Wires) CAF_HitClust.wires.push_back(Wire);
            for(auto Chan : LAr_HitClust.Chans) CAF_HitClust.chans.push_back(Chan);
            for(auto G4ID : LAr_HitClust.G4IDs) CAF_HitClust.G4IDs.push_back(G4ID);
    }
}
