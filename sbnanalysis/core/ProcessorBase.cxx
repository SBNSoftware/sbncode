#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TParameter.h>
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Provenance/SubRunAuxiliary.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "Event.hh"
#include "SubRun.hh"
#include "Loader.hh"
#include "util/Interaction.hh"
#include "ProcessorBase.hh"
#include "ProviderManager.hh"

#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

namespace core {

ProcessorBase::ProcessorBase()
    : fEventIndex(0), fOutputFilename("output.root"), fProviderManager(NULL) {}


ProcessorBase::~ProcessorBase() {}


void ProcessorBase::FillTree() {
  fTree->Fill();
  fEventIndex++;
}


void ProcessorBase::FillRecoTree() {
  fRecoTree->Fill();
}


void ProcessorBase::EventCleanup() {
  fRecoEvent->metadata.Init();
  fEvent->metadata.Init();
  fEvent->truth.clear();
  fEvent->reco.clear();
}


void ProcessorBase::Initialize(char* config) {
  fhicl::ParameterSet* cfg = LoadConfig(config);
  Initialize(cfg);
}


void ProcessorBase::Setup(char* config) {
  fhicl::ParameterSet* cfg = LoadConfig(config);
  Setup(cfg);
}


void ProcessorBase::SetupServices(gallery::Event& ev) {
  if (fProviderManager != NULL) {
    // reset the channels of the back tracker
    if (fProviderManager->GetBackTrackerProvider() != NULL) {
      fProviderManager->GetBackTrackerProvider()->ClearEvent();
      fProviderManager->GetBackTrackerProvider()->PrepSimChannels(ev);
    }

    // reset information in particle inventory
    if (fProviderManager->GetParticleInventoryProvider() != NULL) {
      fProviderManager->GetParticleInventoryProvider()->ClearEvent();
      fProviderManager->GetParticleInventoryProvider()->PrepParticleList(ev);
      fProviderManager->GetParticleInventoryProvider()->PrepMCTruthList(ev);
      fProviderManager->GetParticleInventoryProvider()->PrepTrackIdToMCTruthIndex(ev);
    }
  }
}


void ProcessorBase::Setup(fhicl::ParameterSet* config) {
  // Load configuration file
  if (config) {
    fExperimentID = \
      static_cast<Experiment>(config->get<int>("ExperimentID", kExpOther));
    fFluxTag = { config->get<std::string>("MCFluxTag", "generator") };
    fMCTrackTag = { config->get<std::string>("MCTrackTag", "mcreco") };
    fMCShowerTag = { config->get<std::string>("MCShowerTag", "mcreco") };
    fMCParticleTag = { config->get<std::string>("MCParticleTag", "largeant") };
    fOutputFilename = config->get<std::string>("OutputFile", "output.root");
    fProviderConfig = config->get<std::string>("ProviderConfigFile", "");

    //Get the truth tags 
    fTruthTags = {};
    if (config->has_key("MCTruthTags")) {
      if (config->is_key_to_atom("MCTruthTags")) {
        std::string truth_tag = config->get<std::string>("MCTruthTags");
        fTruthTags = { { truth_tag } };
      }
      else if (config->is_key_to_sequence("MCTruthTags")) {
        std::vector<std::string> truth_tags = \
          config->get<std::vector<std::string> >("MCTruthTags");
        for (size_t i=0; i<truth_tags.size(); i++) {
          fTruthTags.emplace_back(truth_tags[i]);
        }
      }
    }
    else {
      fTruthTags = {"generator"};
    }

    fWriteTree = config->get<bool>("WriteTree", true);
    fWriteRecoTree = config->get<bool>("WriteRecoTree", true);

    // Get the event weight tags (can supply multiple producers)
    fWeightTags = {};
    if (config->has_key("MCWeightTags")) {
      if (config->is_key_to_atom("MCWeightTags")) {
        std::string weight_tag = config->get<std::string>("MCWeightTags");
        fWeightTags = { { weight_tag } };
      }
      else if (config->is_key_to_sequence("MCWeightTags")) {
        std::vector<std::string> weight_tags = \
          config->get<std::vector<std::string> >("MCWeightTags");
        for (size_t i=0; i<weight_tags.size(); i++) {
          fWeightTags.emplace_back(weight_tags[i]);
        }
      }
    }
  }
  // Default -- no config file provided
  else {
    fExperimentID = kExpOther;
    fTruthTags = { "generator", "corskia" };
    fFluxTag = { "generator" };
    fWeightTags = {};
    fMCTrackTag = { "mcreco" };
    fMCShowerTag = { "mcreco" };
    fMCParticleTag = { "largeant" };
    fOutputFilename = "output.root";
  }

  // Set up the provider manager for known experiments
  std::vector<Experiment> exps = ProviderManager::GetValidExperiments();

  if (std::find(exps.begin(), exps.end(), fExperimentID) != exps.end()) {
    std::cout << "fProviderConfig: " << fProviderConfig << std::endl;
    fProviderManager = new ProviderManager(fExperimentID, fProviderConfig);

    // setup the volumes for keeping track of length
    for (auto const &cryo: fProviderManager->GetGeometryProvider()->IterateCryostats()) {
      geo::GeometryCore::TPC_iterator iTPC = fProviderManager->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                      tend = fProviderManager->GetGeometryProvider()->end_TPC(cryo.ID());
      
      // make each cryostat volume a box enclosing all tpc volumes
      double XMin = std::min_element(iTPC, tend, [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
      double YMin = std::min_element(iTPC, tend, [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
      double ZMin = std::min_element(iTPC, tend, [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

      double XMax = std::max_element(iTPC, tend, [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
      double YMax = std::max_element(iTPC, tend, [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
      double ZMax = std::max_element(iTPC, tend, [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

      fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
    }

  }
  else {
    std::cout << "ProcessorBase::Setup: "
              << "Unknown experiment, no ProviderManager is available."
              << std::endl;
  }

  // Open the output file and create the standard event trees
  fOutputFile = TFile::Open(fOutputFilename.c_str(), "recreate");

  fTree = new TTree("sbnana", "SBN Analysis Tree");
  fEvent = new event::Event();
  fTree->Branch("events", &fEvent);
  fReco = &fEvent->reco;

  fRecoTree = new TTree("sbnreco", "SBN Reco Event Tree");
  fRecoEvent = new event::RecoEvent();
  fRecoTree->Branch("reco_events", &fRecoEvent);

  // Create the output subrun tree
  fSubRunTree = new TTree("sbnsubrun", "SBN Analysis Subrun Tree");
  fSubRunTree->AutoSave("overwrite");
  fSubRun = new SubRun();
  fSubRunTree->Branch("subruns", &fSubRun);

  // save the experiment ID
  fExperimentParameter = new TParameter<int>("experiment", fExperimentID);
  fExperimentParameter->Write();

  fTFileName = " ";

}


void ProcessorBase::UpdateSubRuns(gallery::Event& ev) {
  // FIXME: This should use official gallery subrun access once available.
  // N.B. Implementation is fragile and depends on the naming of the subrun
  // producer (generator__GenieGen), can be made a fcl parameter.
  TTree* srtree = (TTree*) ev.getTFile()->Get("SubRuns");

  art::SubRunAuxiliary* sraux = new art::SubRunAuxiliary;
  srtree->SetBranchAddress("SubRunAuxiliary", &sraux);

  for (long i=0; i<srtree->GetEntries(); i++) {
    srtree->GetEntry(i);
    int runid = sraux->run();
    int subrunid = sraux->subRun();
    std::pair<int, int> id = { runid, subrunid };

    // Add subrun if not in cache
    if (fSubRunCache.find(id) == fSubRunCache.end() || fTFileName != ev.getTFile()->GetName()) {
      TLeaf* potLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.totpot");
      float pot = potLeaf ? potLeaf->GetValue() : -1;
      TLeaf* goodpotLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.totgoodpot");
      float goodpot = goodpotLeaf ? goodpotLeaf->GetValue() : -1;
      TLeaf* spillsLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.totspills");
      int spills = spillsLeaf ? spillsLeaf->GetValue() : -1;
      TLeaf* goodspillsLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.goodspills");
      int goodspills = goodspillsLeaf ? goodspillsLeaf->GetValue() : -1;

      if(!potLeaf){
	potLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenGenie.obj.totpot");
	pot = potLeaf ? potLeaf->GetValue() : -1;
      }

      if(!goodpotLeaf){
	 goodpotLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenGenie.obj.totgoodpot");
	 goodpot = goodpotLeaf ? goodpotLeaf->GetValue() : -1;
      }
      if(!spillsLeaf){
	spillsLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenGenie.obj.totspills");
	spills = spillsLeaf ? spillsLeaf->GetValue() : -1;
      }
      if(!goodspillsLeaf){
	goodspillsLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenGenie.obj.goodspills");
	goodspills = goodspillsLeaf ? goodspillsLeaf->GetValue() : -1;
      }
    
      *fSubRun = { runid, subrunid, pot, goodpot, spills, goodspills };
      fSubRunTree->Fill();

      

      fSubRunCache.insert(id);

      std::cout << "Subrun " << runid << "/" << subrunid << " added "
                << "(good POT = " << goodpot << ")"
                << std::endl;
    }
  }
  //Save current file name 
  fTFileName = ev.getTFile()->GetName();
}


void ProcessorBase::Teardown() {
  // Write the standard tree and close the output file
  fOutputFile->cd();
  if (fWriteTree) {
    fTree->Write("sbnana", TObject::kOverwrite);
  }
  if (fWriteRecoTree) {
    fRecoTree->Write("sbnreco", TObject::kOverwrite);
  }
  fSubRunTree->Write("sbnsubrun", TObject::kOverwrite);
  fOutputFile->Purge();
  fOutputFile->Close();
}


void ProcessorBase::BuildEventTree(gallery::Event& ev) {

  // Add any new subruns to the subrun tree
  UpdateSubRuns(ev);

  
  for(int fTT=0; fTT<fTruthTags.size(); ++fTT){
    //    continue;
    // Get MCTruth information
    //auto const& mctruths =						\
    //*ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTags[fTruthTag]);
    
    art::InputTag fTruthTag = fTruthTags[fTT];

    // Get MC truth information
    auto const& mctruths_handle = \
      ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
    const std::vector<simb::MCTruth> &mctruths = *mctruths_handle;
    
    gallery::Handle<std::vector<simb::MCParticle> > mcparticle_list;
    ev.getByLabel(fMCParticleTag, mcparticle_list);
    bool mcparticles_is_valid = mcparticle_list.isValid();

    gallery::Handle<std::vector<simb::GTruth> > gtruths_handle;
    ev.getByLabel(fTruthTag, gtruths_handle);
    bool genie_truth_is_valid = gtruths_handle.isValid();
    
    // get associations from MCTruth information to truth 
    art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truth_to_particles(mctruths_handle, ev, fMCParticleTag);

    // Get MCFlux information
    //    auto const& mcfluxes =					\
    //  *ev.getValidHandle<std::vector<simb::MCFlux> >(fTruthTag);
    //assert(mctruths.size() == mcfluxes.size());

    // Get MCEventWeight information
    std::vector<gallery::Handle<std::vector<evwgh::MCEventWeight> > > wghs;

    if (!fWeightTags.empty()) {
      for (auto const &weightTag: fWeightTags) {
	gallery::Handle<std::vector<evwgh::MCEventWeight> > this_wgh;
	bool hasWeights = ev.getByLabel(weightTag, this_wgh);
	// coherence check
	if (hasWeights) {
	  if(this_wgh->size() == mctruths.size()){wghs.push_back(this_wgh);}
	  //assert(this_wgh->size() == mctruths.size());
	}
	// store the weights
	//	wghs.push_back(this_wgh);
      }
    }
  
    // Get MCFlux information
    gallery::Handle<std::vector<simb::MCFlux> > mcflux_handle;
    ev.getByLabel(fFluxTag, mcflux_handle);

    fTree->GetEntry(fEventIndex);

    // Populate event tree
    fEvent->experiment = fExperimentID;

    auto const& evaux = ev.eventAuxiliary();
    fEvent->metadata.run = evaux.run();
    fEvent->metadata.subrun = evaux.subRun();
    fEvent->metadata.eventID = evaux.event();

    for (size_t i=0; i<mctruths.size(); i++) {
      event::Interaction interaction;
      interaction.index = i;

      auto const& mctruth = mctruths.at(i);
      auto const& mcflux =  mcflux_handle->at(i);

      // TODO: What to do with cosmic MC?
      // For now, ignore them
      //  if (!mctruth.NeutrinoSet()) continue;

      // Combine Weights
      if (!wghs.empty()) {
	size_t wgh_idx = 0;
	// Loop through weight generators, which have a list of weights per truth
	//
	// NOTE: The code allows for multiple different weight generators to produce weights with the same name.
	//       What will happen is that the same-named weights from differetn producers will be placed in the same
	//       vector in the "weights" map below. The order in which weights from different producers are combined
	//       should be consistent from event to event so that correlations are preserved between events. This 
	//       is ensured in the current implementation by making sure that the different weight tags in "wghs"
	//       are always ordered in the same way.
	for (auto const& wgh : wghs) {
	  for (const std::pair<std::string, std::vector<double> >& this_wgh : wgh->at(i).fWeight) {
	    // If we haven't seen this name before, make a new vector with the name
	    if (interaction.weights.count(this_wgh.first) == 0) {
	      interaction.weights.insert({
		  this_wgh.first,
		    std::vector<float>(this_wgh.second.begin(), this_wgh.second.end())
		    });
	    }
	    // If we have seen the name before, append this instance to the last one
	    else {
	      interaction.weights.at(this_wgh.first).insert(
							    interaction.weights.at(this_wgh.first).end(),
							    this_wgh.second.begin(),
							    this_wgh.second.end());
	    }
	  }
	}
      }
  
      if (mcflux_handle.isValid()) {
	const simb::MCFlux& flux = mcflux_handle->at(i);
	interaction.neutrino.parentPDG = flux.fptype;
	interaction.neutrino.parentDecayMode = flux.fndecay;
	interaction.neutrino.parentDecayVtx = \
	  TVector3(flux.fvx, flux.fvy, flux.fvz);
	interaction.neutrino.baseline = flux.fdk2gen + flux.fgen2vtx;
      }

      TLorentzVector q_labframe;

      // get the list of MCParticles to be considered for this interaction
      std::vector<art::Ptr<simb::MCParticle>> this_mcparticle_list; 
      std::vector<const sim::GeneratedParticleInfo *> this_mcparticle_assns;
      if (mcparticles_is_valid) {
	this_mcparticle_list = truth_to_particles.at(i);
	this_mcparticle_assns = truth_to_particles.data(i);
      }

      if (mctruth.NeutrinoSet()) {
	// Neutrino
	const simb::MCNeutrino& nu = mctruth.GetNeutrino();
	interaction.neutrino.isnc =   nu.CCNC()  && (nu.Mode() != simb::kWeakMix);
	interaction.neutrino.iscc = (!nu.CCNC()) && (nu.Mode() != simb::kWeakMix);
	interaction.neutrino.pdg = nu.Nu().PdgCode();
	interaction.neutrino.initpdg = mcflux.fntype;
	interaction.neutrino.targetPDG = nu.Target();
	interaction.neutrino.genie_intcode = nu.Mode();
	interaction.neutrino.bjorkenX = nu.X();
	interaction.neutrino.inelasticityY = nu.Y();
	interaction.neutrino.Q2 = nu.QSqr();
	interaction.neutrino.w = nu.W();
	interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();
	interaction.neutrino.momentum = nu.Nu().EndMomentum().Vect();
	interaction.neutrino.position = nu.Nu().Position().Vect();

	// Primary lepton
	const simb::MCParticle& lepton = nu.Lepton();
	interaction.lepton.pdg = lepton.PdgCode();
	interaction.lepton.energy = lepton.Momentum(0).Energy();
	interaction.lepton.momentum = lepton.Momentum(0).Vect();
	interaction.lepton.start = lepton.Position(0).Vect();
	interaction.lepton.status_code = lepton.StatusCode();
	interaction.lepton.is_primary = (lepton.Process() == "primary");

	// match the MCTruth particle to the MCParticle list -- only the list has trajectory information
	const simb::MCParticle* lepton_traj = NULL;
	for (int iparticle = 0; iparticle < this_mcparticle_list.size(); iparticle++) {
	  if (this_mcparticle_assns[iparticle]->hasGeneratedParticleIndex() &&
	      mctruth.GetParticle(this_mcparticle_assns[iparticle]->generatedParticleIndex()).TrackId() == lepton.TrackId()) {

	    lepton_traj = this_mcparticle_list[iparticle].get();
	    break;
	  }
	}
    
	// TODO: why aren't there matches sometimes?
	if (lepton_traj != NULL) {
	  interaction.lepton.length = util::MCParticleLength(*lepton_traj);
	  if (fProviderManager != NULL) {
	    interaction.lepton.contained_length = util::MCParticleContainedLength(*lepton_traj, fActiveVolumes);
	  }
	  else {
	    interaction.lepton.contained_length = event::kUnfilled;
	  }
	  // also get the end point from the trajectory
	  interaction.lepton.end = lepton_traj->EndPosition().Vect();
	}
	else {
	  interaction.lepton.contained_length = 0;
	  interaction.lepton.length = 0;
	  // if we couldn't find a trajectory, set end to start
	  interaction.lepton.end = interaction.lepton.start;
	}

	q_labframe = nu.Nu().EndMomentum() - lepton.Momentum(0);
	interaction.neutrino.q0_lab = q_labframe.E();
	interaction.neutrino.modq_lab = q_labframe.P();
      }

      // Get CCQE energy from lepton info
      interaction.neutrino.eccqe = \
	util::ECCQE(interaction.lepton.momentum, interaction.lepton.energy);

      // Hadronic system
      for (int iparticle=0; iparticle<this_mcparticle_list.size(); iparticle++) {
	event::FinalStateParticle fsp;
	const simb::MCParticle& particle = *this_mcparticle_list[iparticle];

	if (particle.Process() != "primary") {
	  continue;
	}

	fsp.pdg = particle.PdgCode();
	fsp.energy = particle.Momentum(0).Energy();
	fsp.momentum = particle.Momentum(0).Vect();
	fsp.start = particle.Position(0).Vect();
	fsp.status_code = particle.StatusCode();
	fsp.is_primary = (particle.Process() == "primary");
	fsp.length = util::MCParticleLength(particle);
	if (fProviderManager != NULL) {
	  fsp.contained_length = util::MCParticleContainedLength(particle, fActiveVolumes);
	}
	else {
	  fsp.contained_length = event::kUnfilled;
	}
	// also get the end point from the trajectory
	fsp.end = particle.EndPosition().Vect();

	interaction.finalstate.push_back(fsp);
      }

      interaction.nfinalstate = interaction.finalstate.size();

      // GENIE specific
      if (genie_truth_is_valid) {
	auto const& gtruth = gtruths_handle->at(i);
	TLorentzVector q_nucframe(q_labframe);
	// This nucleon momentum should be added to MCNeutrino so we don't
	// have to rely on GTruth
	const TLorentzVector& nucP4 = gtruth.fHitNucP4;
	TVector3 nuc_boost(nucP4.BoostVector());
	q_nucframe.Boost(nuc_boost);
	interaction.neutrino.modq = q_nucframe.P();
	interaction.neutrino.q0 = q_nucframe.E();
      }
      fEvent->truth.push_back(interaction);
    }
  }

  fEvent->ntruth = fEvent->truth.size();
}

}  // namespace core

