#include "SetEvent.h"
#include "PIDDataMaker.h"

#include "../NumuReco/TrackAlgo.h"

void ana::SBNOsc::PIDDataMaker::Initialize(fhicl::ParameterSet *config) {
  fCuts.Initialize(config->get<fhicl::ParameterSet>("Cuts"), fProviderManager->GetGeometryProvider());
  fOutputFile = new TFile(config->get<std::string>("OutputFile", "output.root").c_str(), "CREATE");

  fPIDData = new numu::flat::PIDData;
  fPIDTree = new TTree("pid", "pid"); 
  fPIDTree->Branch("pid", &fPIDData);
  fRecoEvent = NULL;
}

void ana::SBNOsc::PIDDataMaker::FileSetup(TFile *f, TTree *eventTree) {
  eventTree->SetBranchAddress("reco_event", &fRecoEvent);
}

void ana::SBNOsc::PIDDataMaker::ProcessEvent(const event::Event *ev) {
  unsigned index = WorkerID();
  ana::SBNOsc::SetEvent(*fRecoEvent, *ev, fCuts, numu::fOverlay);
  for (const numu::RecoInteraction &interaction: fRecoEvent->reco) {
    int truth_index = interaction.slice.truth.interaction_id;
    if (truth_index < 0 /* cosmic */) continue;
    if (!fCuts.InFV(interaction.position) /* non-fiducial */) continue;

    for (unsigned ID: interaction.PrimaryTracks(fRecoEvent->tracks)) {
      const numu::RecoTrack &track = fRecoEvent->tracks.at(ID);
      if (!fRecoEvent->particles.count(track.truth.GetPrimaryMatchID()) /* no valid truth match */) continue;
      const numu::TrueParticle &particle = fRecoEvent->particles.at(track.truth.GetPrimaryMatchID());
      
      fPIDData->true_pdg = abs(particle.pdgid); 
      fPIDData->trueP = particle.start_momentum.Mag();
      fPIDData->purity = track.truth.Purity();
      fPIDData->completion = track.truth.matches[0].energy / particle.deposited_energy;
      fPIDData->length = track.length;
      fPIDData->phi = track.phi;
      fPIDData->theta = track.theta;
      fPIDData->chi2_muon = track.chi2_muon;
      fPIDData->chi2_proton = track.chi2_proton;
      fPIDData->chi2_pion = track.chi2_pion;
      fPIDData->crt_hit_distance = -1;
      if (track.crt_match.hit_match.present) fPIDData->crt_hit_distance = track.crt_match.hit_match.distance;
      fPIDData->contained = track.is_contained;
      fPIDData->n_trk_daughters = 0;
      fPIDData->n_shr_daughters = 0;
      for (unsigned id: interaction.slice.particles.at(ID).daughters) {
        if (fRecoEvent->tracks.count(id)) fPIDData->n_trk_daughters ++;
        else fPIDData->n_shr_daughters ++;
      }

      fPIDTree->Fill();
    }
  }
}

void ana::SBNOsc::PIDDataMaker::Finalize() {
  fOutputFile->cd();
  fPIDTree->Write();
}



DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::PIDDataMaker);
