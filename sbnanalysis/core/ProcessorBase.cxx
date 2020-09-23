#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TParameter.h>
#include <unordered_map>
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
#include "fhiclcpp/ParameterSetID.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "cetlib/sqlite/query_result.h"
#include "cetlib/sqlite/select.h"

#include "canvas/Persistency/Provenance/ParameterSetBlob.h"
#include "canvas/Persistency/Provenance/ParameterSetMap.h"
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "art_root_io/RootDB/tkeyvfs.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "Event.hh"
#include "SubRun.hh"
#include "Loader.hh"
#include "util/Interaction.hh"
#include "ProcessorBase.hh"
#include "ProviderManager.hh"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
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


void ProcessorBase::Setup(fhicl::ParameterSet* config) {
  // Load configuration file
  if (config) {
    fExperimentID = \
      static_cast<Experiment>(config->get<int>("ExperimentID", kExpOther));
    fTruthTag = { config->get<std::string>("MCTruthTag", "generator") };
    fGeneratorProcess = config->get<std::string>("GeneratorProcess", "GenieGen");
    fFluxTag = { config->get<std::string>("MCFluxTag", "generator") };
    fMCTrackTag = { config->get<std::string>("MCTrackTag", "mcreco") };
    fMCShowerTag = { config->get<std::string>("MCShowerTag", "mcreco") };
    fMCParticleTag = { config->get<std::string>("MCParticleTag", "largeant") };
    fOutputFilename = config->get<std::string>("OutputFile", "output.root");
    fProviderConfig = config->get<std::string>("ProviderConfigFile", "");

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
    fTruthTag = { "generator" };
    fGeneratorProcess = "GenieGen";
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

  // Create the file metadata tree
  fFileMetaTree = new TTree("sbnfilemeta", "SBN Analysis File Metadata Tree");
  fFileMeta = new FileMeta();
  fFileMetaTree->Branch("filemeta", &fFileMeta);

  // save the experiment ID
  fExperimentParameter = new TParameter<int>("experiment", fExperimentID);
  fExperimentParameter->Write();

  // setup the connection to the sqlite database
  tkeyvfs_init();
}

void ProcessorBase::UpdateSubRuns(gallery::Event& ev) {
  // FIXME: HACK
  static TFile *last_tfile = NULL; // only set at initialization

  TFile *this_file = ev.getTFile();
  // it's a new file!
  if (this_file != last_tfile) {
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

      std::string totpot_str = "sumdata::POTSummary_generator__" + fGeneratorProcess + ".obj.totpot";
      std::string totgoodpot_str = "sumdata::POTSummary_generator__" + fGeneratorProcess + ".obj.totgoodpot";
      std::string totspills_str = "sumdata::POTSummary_generator__" + fGeneratorProcess + ".obj.totspills";
      std::string goodspills_str = "sumdata::POTSummary_generator__" + fGeneratorProcess + ".obj.goodspills";

      TLeaf* potLeaf = srtree->GetLeaf(totpot_str.c_str());
      float pot = potLeaf ? potLeaf->GetValue() : -1;
      TLeaf* goodpotLeaf = srtree->GetLeaf(totgoodpot_str.c_str());
      float goodpot = goodpotLeaf ? goodpotLeaf->GetValue() : -1;
      TLeaf* spillsLeaf = srtree->GetLeaf(totspills_str.c_str());
      int spills = spillsLeaf ? spillsLeaf->GetValue() : -1;
      TLeaf* goodspillsLeaf = srtree->GetLeaf(goodspills_str.c_str());
      int goodspills = goodspillsLeaf ? goodspillsLeaf->GetValue() : -1;

      std::cout << "Subrun " << runid << "/" << subrunid << " added "
                << "(good POT = " << goodpot << ")"
                << std::endl;
      *fSubRun = { runid, subrunid, pot, goodpot, spills, goodspills };
      fSubRunTree->Fill();
    }
  }
  // set the file for next time
  last_tfile = ev.getTFile();
}

void ProcessorBase::UpdateFileMeta(gallery::Event& ev) {
  // FIXME: this is also a bit of a hack
  static TFile *last_tfile = NULL; // only set at initialization

  TFile *this_file = ev.getTFile();
  // it's a new file!
  if (this_file != last_tfile) {
    unsigned n_gen_event = 0;
    unsigned n_event = 0;

    // get the SQLite thingy
    art::SQLite3Wrapper sqliteDB(this_file, "RootFileDB", SQLITE_OPEN_READONLY);
    sqlite3 *db = sqliteDB;
    // add to the registry and parse yourself
    fhicl::ParameterSetRegistry::importFrom(db);
    fhicl::ParameterSetRegistry::stageIn();

    // TODO: make this better
    // 
    std::unordered_map<fhicl::ParameterSetID, fhicl::ParameterSet, fhicl::detail::HashParameterSetID> registry;
    using namespace cet::sqlite;
    query_result<std::string, std::string> entriesToStageIn;
    entriesToStageIn << select("*").from(db, "ParameterSets");

    cet::transform_all(entriesToStageIn,
                     std::inserter(registry, std::begin(registry)),
                     [](auto const& row) {
                       auto const& [idString, psBlob] = row;
                       fhicl::ParameterSet pset;
                       fhicl::make_ParameterSet(psBlob, pset);
                       return std::make_pair(fhicl::ParameterSetID{idString}, pset);
                     });

    for (auto const& pr: registry) {
      const fhicl::ParameterSet pset = pr.second;
      if (pset.has_key("source") && pset.has_key("source.maxEvents") && pset.has_key("source.module_type")) {
        int max_events = pset.get<int>("source.maxEvents");
        std::string module_type = pset.get<std::string>("source.module_type");
        // this will show up once for every input file
        if (module_type == "EmptyEvent") {
          n_gen_event += max_events;
          std::cout << "FILEMETA: " << max_events << std::endl;
        }
      }
    }
    n_event = ev.numberOfEventsInFile();
    std::cout << "FILEMETA nevent: " << n_event << std::endl;
    std::cout << "FILEMETA ngenevent: " << n_gen_event << std::endl;
    *fFileMeta = {n_event, n_gen_event};
    fFileMetaTree->Fill();
  }
  // set the file for next time
  last_tfile = ev.getTFile();
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
  fFileMetaTree->Write("sbnfilemeta", TObject::kOverwrite);
  fOutputFile->Purge();
  fOutputFile->Close();
}

void ProcessorBase::SetupServices(gallery::Event& ev) {
  if (fProviderManager != NULL) {
    fProviderManager->SetupServices(ev);
  }
}

const simb::MCParticle *Genie2G4MCParticle(
  const simb::MCParticle &genie_part, 
  const simb::MCTruth &mctruth,
  const std::vector<art::Ptr<simb::MCParticle>> &g4_mcparticles, 
  const std::vector<const sim::GeneratedParticleInfo *> infos) {

  const simb::MCParticle *ret = NULL;
  for (int iparticle = 0; iparticle < g4_mcparticles.size(); iparticle++) {
    if (infos[iparticle]->hasGeneratedParticleIndex() &&
        infos[iparticle]->generatedParticleIndex() < mctruth.NParticles() && // TODO: why is this number sometimes bigger than the number of particles?
        mctruth.GetParticle(infos[iparticle]->generatedParticleIndex()).TrackId() == genie_part.TrackId() &&
        g4_mcparticles[iparticle]->Process() == "primary" /* TODO: will have to remove this restriction to include secondary particles*/) {

      // if a genie particle re-scatters in g4 and makes more particles, then multiple g4 particles can match to a 
      // genie particle. Thus, we also check that the start location of the associated genie particle matches the g4
      // and that the pdgid matches (to be on the safe side)
      //
      // Note that this should be accounted for by requiring the process to be primary. This is a bit of redundancy.
      const simb::MCParticle& matched_genie_particle = mctruth.GetParticle(infos[iparticle]->generatedParticleIndex());
      if ((matched_genie_particle.Position().Vect() - g4_mcparticles[iparticle]->Position().Vect()).Mag() < 1e-4 &&
        matched_genie_particle.PdgCode() == g4_mcparticles[iparticle]->PdgCode()) {

        // this should only be true for one particle
        assert(ret == NULL);
        ret = g4_mcparticles[iparticle].get();
      }
    }
  }
  return ret;
}

void ProcessorBase::BuildEventTree(gallery::Event& ev) {
  // Add any new subruns to the subrun tree
  UpdateSubRuns(ev);
  UpdateFileMeta(ev);

  // Get MC truth information
  gallery::Handle<std::vector<simb::MCTruth>> mctruths_handle;
  bool mctruth_is_valid = ev.getByLabel(fTruthTag, mctruths_handle);

  gallery::Handle<std::vector<simb::MCParticle> > mcparticle_list;
  ev.getByLabel(fMCParticleTag, mcparticle_list);
  bool mcparticles_is_valid = mcparticle_list.isValid();

  gallery::Handle<std::vector<simb::GTruth> > gtruths_handle;
  ev.getByLabel(fTruthTag, gtruths_handle);
  bool genie_truth_is_valid = gtruths_handle.isValid();

  // Get MCEventWeight information
  std::vector<gallery::Handle<std::vector<evwgh::MCEventWeight> > > wghs;

  if (!fWeightTags.empty()) {
    for (auto const &weightTag: fWeightTags) {
      gallery::Handle<std::vector<evwgh::MCEventWeight> > this_wgh;
      bool hasWeights = ev.getByLabel(weightTag, this_wgh);
      // coherence check
      if (hasWeights && mctruth_is_valid) {
        assert(this_wgh->size() == mctruths_handle->size());
      }
      // store the weights
      wghs.push_back(this_wgh);
    }
  }

  // Get MCFlux information
  gallery::Handle<std::vector<simb::MCFlux> > mcflux_handle;
  bool flux_is_valid = ev.getByLabel(fFluxTag, mcflux_handle);
  // Check MCFlux information
  if (flux_is_valid && mctruth_is_valid) {
    assert(mctruths_handle->size() == mcflux_handle->size());
  }


  fTree->GetEntry(fEventIndex);

  // Populate event tree
  fEvent->experiment = fExperimentID;

  auto const& evaux = ev.eventAuxiliary();
  fEvent->metadata.run = evaux.run();
  fEvent->metadata.subrun = evaux.subRun();
  fEvent->metadata.eventID = evaux.event();
  fEvent->metadata.fileEntry = ev.fileEntry();

  if (mctruth_is_valid) {
    // get associations from MCTruth information to truth 
    art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truth_to_particles(mctruths_handle, ev, fMCParticleTag);
    unsigned n_truth = mctruths_handle->size();
    for (size_t i=0; i<n_truth; i++) {

      event::Interaction interaction;
      interaction.index = i;

      auto const& mctruth = mctruths_handle->at(i);

      // TODO: What to do with cosmic MC?
      // For now, ignore them
      if (!mctruth.NeutrinoSet()) continue;

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
        interaction.neutrino.initpdg = flux.fntype;
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
	const simb::MCParticle* lepton_traj = Genie2G4MCParticle(lepton, mctruth, this_mcparticle_list, this_mcparticle_assns);

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
      //
      // We need to look at the Genie MCParticles because we may want
      // to look at truth information like (e.g.) intermeidate-state pions
      // which are re-absorbed in the nucleus. 
      //
      // Whenever we can though, try to match the true "gen" information to g4
      // information
      for (int iparticle=0; iparticle<mctruth.NParticles(); iparticle++) {
        const simb::MCParticle& particle = mctruth.GetParticle(iparticle);
        event::FinalStateParticle fsp;
        if (particle.Process() != "primary") {
          continue;
        }
	fsp.pdg = particle.PdgCode();
	fsp.energy = particle.Momentum(0).Energy();
	fsp.momentum = particle.Momentum(0).Vect();
	fsp.start = particle.Position(0).Vect();
	fsp.status_code = particle.StatusCode();
	fsp.is_primary = (particle.Process() == "primary");
	
	fsp.length = 0;
	fsp.contained_length = 0;
	// length information needs G4
	const simb::MCParticle *traj = Genie2G4MCParticle(particle, mctruth, this_mcparticle_list, this_mcparticle_assns);
	if (traj != NULL) {
          fsp.length = util::MCParticleLength(*traj);
          if (fProviderManager != NULL) {
            fsp.contained_length = util::MCParticleContainedLength(*traj, fActiveVolumes);
          }
          else {
            fsp.contained_length = event::kUnfilled;
          }
        }
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

