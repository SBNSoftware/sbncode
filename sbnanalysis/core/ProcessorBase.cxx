#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "gallery/FindMaker.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Provenance/SubRunAuxiliary.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "larsim/MCCheater/ParticleInventoryService.h" 
#include "larsim/MCCheater/BackTrackerService.h" 
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "Event.hh"
#include "SubRun.hh"
#include "Loader.hh"
#include "util/Interaction.hh"
#include "util/RecoUtil.hh"
#include "ProcessorBase.hh"
#include "ProviderManager.hh"

namespace core {

  ProcessorBase::ProcessorBase()
    : fEventIndex(0), fOutputFilename("output.root"), fProviderManager(NULL), fMCSFitter(NULL) {}

  ProcessorBase::~ProcessorBase() {
    delete fMCSFitter;
  }

  void ProcessorBase::FillTree() {
    fTree->Fill();
    fEventIndex++;
  }

  void ProcessorBase::EventCleanup() {
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
    // With configuration file provided
    if (config) {
      fExperimentID = \
                      static_cast<Experiment>(config->get<int>("ExperimentID", kExpOther));
      fTruthTag                = config->get<std::string>("MCTruthTag", "generator");
      fFluxTag                 = config->get<std::string>("MCFluxTag", "generator");
      fMCTrackTag              = config->get<std::string>("MCTrackTag", "mcreco");
      fMCShowerTag             = config->get<std::string>("MCShowerTag", "mcreco");
      fMCParticleTag           = config->get<std::string>("MCParticleTag", "largeant");
      fHitTag                  = config->get<std::string>("HitTag", "gaushit");
      fPFParticleTag           = config->get<std::string>("PFParticleTag", "pandora");
      fRecoTrackTag            = config->get<std::string>("RecoTrackTag", "pandoraTrack");
      fRecoShowerTag           = config->get<std::string>("RecoShowerTag", "emshower");
      fVertexTag               = config->get<std::string>("VertexTag", "pandora");
      fRecoTrackCalorimetryTag = config->get<std::string>("RecoTrackCalorimetryTag", "pandoraCalo");
      fRecoTrackParticleIDTag  = config->get<std::string>("RecoTrackParticleIDTag", "pandoraPid");
      fOutputFilename          = config->get<std::string>("OutputFile", "output.root");
      fProviderConfig          = config->get<std::string>("ProviderConfigFile", "");
      fMCSFitterPID            = config->get<int>("pIdHypothesis");
      fMCSFitterNumSegments    = config->get<int>("minNumSegments");
      fMCSFitterHitsSegment    = config->get<int>("minHitsPerSegment");
      fMCSFitterElossSteps     = config->get<int>("nElossSteps");
      fMCSFitterElossMode      = config->get<int>("eLossMode");
      fMCSFitterLengthSegment  = config->get<double>("segmentLength");
      fMCSFitterPMin           = config->get<double>("pMin");
      fMCSFitterPMax           = config->get<double>("pMax");
      fMCSFitterPStep          = config->get<double>("pStep");
      fMCSFitterAngResol       = config->get<double>("angResol");
      fMCSFitter = new trkf::TrajectoryMCSFitter(fMCSFitterPID, 
                                                 fMCSFitterNumSegments,
                                                 fMCSFitterLengthSegment,
                                                 fMCSFitterHitsSegment,
                                                 fMCSFitterElossSteps, 
                                                 fMCSFitterElossMode,  
                                                 fMCSFitterPMin,   
                                                 fMCSFitterPMax,        
                                                 fMCSFitterPStep,       
                                                 fMCSFitterAngResol);


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
      fFluxTag = { "generator" };
      fWeightTags = {};
      fMCTrackTag = { "mcreco" };
      fMCShowerTag = { "mcreco" };
      fMCParticleTag = { "largeant" };
      fHitTag        = { "gaushit" };
      fPFParticleTag = { "pandora" };
      fRecoTrackTag = { "pandoraTrack" };
      fRecoShowerTag = { "emshower" };
      fVertexTag     = { "pandora" };
      fRecoTrackCalorimetryTag = { "pandoraCalo" };
      fRecoTrackParticleIDTag = { "pandoraPid" };
      fOutputFilename = "output.root";
    }

    // Set up the provider manager for known experiments
    std::vector<Experiment> exps = ProviderManager::GetValidExperiments();
    if (std::find(exps.begin(), exps.end(), fExperimentID) != exps.end()) {
      fProviderManager = new ProviderManager(fExperimentID, fProviderConfig);
    }
    else {
      std::cout << "ProcessorBase::Setup: "
        << "Unknown experiment, no ProviderManager is available."
        << std::endl;
    }

    // Open the output file and create the standard event tree
    fOutputFile = TFile::Open(fOutputFilename.c_str(), "recreate");

    fTree = new TTree("sbnana", "SBN Analysis Tree");
    fTree->AutoSave("overwrite");
    fEvent = new Event();
    fTree->Branch("events", &fEvent);
    fReco = &fEvent->reco;

    // Create the output subrun tree
    fSubRunTree = new TTree("sbnsubrun", "SBN Analysis Subrun Tree");
    fSubRunTree->AutoSave("overwrite");
    fSubRun = new SubRun();
    fSubRunTree->Branch("subruns", &fSubRun);
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
      if (fSubRunCache.find(id) == fSubRunCache.end()) {
        TLeaf* potLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.totpot");
        double pot = potLeaf ? potLeaf->GetValue() : -1;
        TLeaf* goodpotLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.totgoodpot");
        double goodpot = goodpotLeaf ? goodpotLeaf->GetValue() : -1;
        TLeaf* spillsLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.totspills");
        int spills = spillsLeaf ? spillsLeaf->GetValue() : -1;
        TLeaf* goodspillsLeaf = srtree->GetLeaf("sumdata::POTSummary_generator__GenieGen.obj.goodspills");
        int goodspills = goodspillsLeaf ? goodspillsLeaf->GetValue() : -1;

        *fSubRun = { runid, subrunid, pot, goodpot, spills, goodspills };
        fSubRunTree->Fill();

        fSubRunCache.insert(id);

        std::cout << "Subrun " << runid << "/" << subrunid << " added "
          << "(good POT = " << goodpot << ")"
          << std::endl;
      }
    }
  } // UpdateSubRuns

  void ProcessorBase::Teardown() {
    // Write the standard tree and close the output file
    fOutputFile->cd();
    fTree->Write("sbnana", TObject::kOverwrite);
    fSubRunTree->Write("sbnsubrun", TObject::kOverwrite);
    fOutputFile->Close();
  }// Teardown
  
  void ProcessorBase::SetupServices(gallery::Event& ev) {
    if (fProviderManager != NULL) {
      /*
      // reset the channels of the back tracker
      fProviderManager->GetBackTrackerProvider()->ClearEvent();
      fProviderManager->GetBackTrackerProvider()->PrepSimChannels(ev);

      // reset information in particle inventory
      fProviderManager->GetParticleInventoryProvider()->ClearEvent();
      fProviderManager->GetParticleInventoryProvider()->PrepParticleList(ev);
      fProviderManager->GetParticleInventoryProvider()->PrepMCTruthList(ev);
      fProviderManager->GetParticleInventoryProvider()->PrepTrackIdToMCTruthIndex(ev);
      */
    }
  } // SetupServices

  void ProcessorBase::BuildEventTree(gallery::Event& ev) {
    // Add any new subruns to the subrun tree
    UpdateSubRuns(ev);

    // Get vector of MCTruth information
    auto const &mctruths    = *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
    int mctruthssize        = mctruths.size();

    // Get vector of MCParticle information
    auto const &mcparticles = *ev.getValidHandle<std::vector<simb::MCParticle> >(fMCParticleTag);
    int mcparticlessize     = mcparticles.size();

    // Get the hit handle
    gallery::Handle< std::vector< recob::Hit > > hits;
    std::vector< art::Ptr< recob::Hit > > hit_vector;
    if(ev.getByLabel(fHitTag, hits)) art::fill_ptr_vector(hit_vector, hits);
    int hitssize            = hits->size();

    // Get the PFParticle handle
    auto const pfparticles  = ev.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleTag);
    int pfparticlessize     = pfparticles->size();

    // Get the Track handle
    auto const recotracks   = ev.getValidHandle<std::vector<recob::Track>>(fRecoTrackTag);
    int recotrackssize      = recotracks->size();

    // Get vector of Shower information
    auto const recoshowers  = ev.getValidHandle<std::vector<recob::Shower> >(fRecoShowerTag);
    int showersize          = recoshowers->size();

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
        if (hasWeights) {
          assert(this_wgh->size() == mctruths.size());
        }
        // store the weights
        wghs.push_back(this_wgh);
      }// Weight tags loop
    }// Weight tags

    // Get MCFlux information
    gallery::Handle<std::vector<simb::MCFlux> > mcflux_handle;
    ev.getByLabel(fFluxTag, mcflux_handle);

    fTree->GetEntry(fEventIndex);

    // Populate event tree
    fEvent->experiment = fExperimentID;

    // Accessing the geometry information needed
    //
    // Variables to get from fProviderManager()->GetGeometryProvider()
    //    Detector borders in x, y and z
    //    Print these
    if(!fProviderManager){
      std::cerr << " Error: No provider manager available " << std::endl;
      exit(1);
    }
    const geo::GeometryCore* geom = fProviderManager->GetGeometryProvider();

    std::vector<double> minx, miny, minz, maxx, maxy, maxz;
    minx.clear();
    miny.clear();
    minz.clear();
    maxx.clear();
    maxy.clear();
    maxz.clear();

    geo::GeometryCore::TPC_id_iterator iTPC = geom->begin_TPC_id(),
      tend = geom->end_TPC_id();
    while (iTPC != tend) {
      // the TPC descriptor object
      const geo::TPCGeo* tpcgeo = iTPC.get();
      // the cryostat the TPC is in
      //geo::CryostatGeo const& Cryo = geom->Cryostat(*iTPC);

      minx.push_back(tpcgeo->MinX()); maxx.push_back(tpcgeo->MaxX());
      miny.push_back(tpcgeo->MinY()); maxy.push_back(tpcgeo->MaxY());
      minz.push_back(tpcgeo->MinZ()); maxz.push_back(tpcgeo->MaxZ());
      ++iTPC;
    } // Loop over the TPCs
    // Populate event tree
    if(!mctruthssize){
      std::cerr << " Error: No MCTruth objects available " << std::endl;
      return;
    }
    // Vectors to hold information about the 'primaries' to determine whether they are
    // neutrinos or cosmics and how many were reconstructed
    std::vector<int> neutrinos;
    std::vector<int>::iterator it;

    // Loop over PFParticles, find neutrinos and add their IDs to the list above
    // if they are not already listed
    for(int i = 0; i < pfparticlessize; ++i){

      // Get the current PFParticle
      auto const pfparticle = pfparticles->at(i);

      // Make sure we are looking at primary final state particles from the neutrino
      if(pfparticle.IsPrimary()){
        // Check if we have reconstructed a neutrino, if not, quit the event
        // This will need to change eventually
        if(pfparticle.PdgCode() != 12 && pfparticle.PdgCode() != 14) return;
        it = std::find(neutrinos.begin(), neutrinos.end(),  pfparticle.Self());
        if(it != neutrinos.end()) continue;
        neutrinos.push_back(pfparticle.Self());
      } // Neutrino check
    } // pfparticles

    // If the number of reconstructed neutrinos does not match the true number
    // of neutrinos, quit. This will get messy as the reconstruction has failed
    if(neutrinos.size() != mctruthssize) {
      std::cerr << " Number of truth objects, " << mctruthssize << " doesn't match the number of reconstructed neutrino objects, " << neutrinos.size() << ": Skipping " << std::endl;
      return;
    }

    for (size_t i=0; i<mctruthssize; i++) {
      Event::Interaction interaction;

      auto const& mctruth = mctruths.at(i);
      //std::vector< art::Ptr<simb::MCParticle> > mcp_assn = fmcp.at(i);

      // TODO: What to do with cosmic MC?
      // For now, ignore them
      if (!mctruth.NeutrinoSet()) continue;

      // Combine Weights
      if (!wghs.empty()) {
        for (auto const &wgh: wghs) {
          // Insert the weights for each individual EventWeight object into the 
          // Event class "master" weight list
          interaction.weights.insert(wgh->at(i).fWeight.begin(), wgh->at(i).fWeight.end());
        }
      }

      if (mcflux_handle.isValid()) {
        const simb::MCFlux& flux = mcflux_handle->at(i);
        interaction.neutrino.parentPDG = flux.fptype;
        interaction.neutrino.parentDecayMode = flux.fndecay;
        interaction.neutrino.parentDecayVtx = \
                                              TVector3(flux.fvx, flux.fvy, flux.fvz);
      }

      TLorentzVector q_labframe;

      if (!mctruth.NeutrinoSet()) {
        std::cerr << " Error: No neutrino information available "<< std::endl;
        return;
      }
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

      q_labframe = nu.Nu().EndMomentum() - lepton.Momentum(0);
      interaction.neutrino.q0_lab = q_labframe.E();
      interaction.neutrino.modq_lab = q_labframe.P();

      // Get CCQE energy from lepton info
      interaction.neutrino.eccqe = \
                                   util::ECCQE(interaction.lepton.momentum, interaction.lepton.energy);

      // Hadronic system
      // Populate map of mctruth IDs and Event::Interactions
      // This will be used to compare PFParticles to MCParticles for the
      // sake of matching MC Interactions to true interactions in an event
      // Also define an iterator
      std::map<int, Event::Interaction> mcinteraction_map;
      std::map<int, Event::Interaction>::const_iterator mcinteraction_map_iterator;

      // Hadronic system
      if(mcparticlessize == 0) continue;
      for(const simb::MCParticle particle : mcparticles){
        Event::FinalStateParticle fsp;
        if (particle.Process() != "primary" || particle.PdgCode() >= 1000018039) continue;

        std::vector< art::Ptr< recob::Hit > > trk_hits = fProviderManager->GetBackTrackerProvider()->TrackIdToHits_Ps(particle.TrackId(), hit_vector);

        fsp.mc_id     = particle.TrackId();
        fsp.pdg       = particle.PdgCode();
        fsp.energy    = particle.Momentum(0).Energy();
        fsp.momentum  = particle.Momentum(0).Vect();
        fsp.vertex    = particle.Position().Vect();
        fsp.end       = particle.EndPosition().Vect();
        fsp.trans_mom = particle.Pt();
        fsp.hits      = trk_hits.size(); 

        interaction.finalstate.push_back(fsp);
      } // mcparticles

      // GENIE specific
      if (genie_truth_is_valid) {
        auto const& gtruth = gtruths_handle->at(i);
        TLorentzVector q_nucframe(q_labframe);
        // This nucleon momentum should be added to MCNeutrino so we don't
        // have to rely on GTruth
        const TLorentzVector& nucP4 = gtruth.fHitNucP4;
        TVector3 nuc_boost(nucP4.BoostVector());
        q_nucframe.Boost(nuc_boost);
        interaction.neutrino.modq     = q_nucframe.P();
        interaction.neutrino.q0       = q_nucframe.E();
        interaction.neutrino.gscatter = gtruth.fGscatter;
        interaction.neutrino.gint     = gtruth.fGint;
      } // GENIE

      mcinteraction_map.insert(std::pair<int, Event::Interaction>(i, interaction));
      fEvent->truth.push_back(interaction);
    }// mctruths

    // Get the track and vertex associations to the PFParticles
    art::FindManyP< recob::Track > fmtrk(pfparticles, ev, fRecoTrackTag);
    art::FindManyP< recob::Vertex > fvtx( pfparticles, ev, fVertexTag);
    // If we eventually use the Pandora Shower objects instead
    //art::FindManyP< recob::Shower > fmshw(pfparticles, ev, fRecoShowerTag);

    // Now loop over the neutrino IDs and the pfparticles 
    // Make a new RecoInteraction for each NeutrinoID 
    // fill the finalstatereconstructedparticles for each of its daughters
    for(int i = 0; i < neutrinos.size(); ++i){

      Event::RecoInteraction reconstructed_interaction;
      std::vector< art::Ptr<recob::Vertex> > vtx_assn = fvtx.at(neutrinos[i]);
      if(vtx_assn.size()  > 1 || vtx_assn.size() == 0) return;

      // Set array to be current vertex position
      recob::Vertex::Point_t vtx = vtx_assn[0]->position();
      TVector3 vertex(vtx.X(), vtx.Y(), vtx.Z());
      reconstructed_interaction.reco_vertex = vertex;

      // Loop over PFParticles, find neutrinos and add their IDs to the list above
      // if they are not already listed
      for(int i = 0; i < pfparticlessize; ++i){
        // Get the current PFParticle
        auto const pfparticle = pfparticles->at(i);
        // If the parent of the current particle is the neutrino

        if(pfparticle.Parent() == neutrinos[i]){
          Event::FinalStateReconstructedParticle fsrp;

          // Let's check if the inheritance is correct by filling a member of 
          // FinalStateParticle
          fsrp.pdg = pfparticle.PdgCode();

          // Get the track associations to get all the other associations
          std::vector< art::Ptr<recob::Track> > trk_assn = fmtrk.at(pfparticle.Self());

          // If the tracks exist, get the hits, particleID and calorimetry associations to them
          if(trk_assn.size()) {

            art::FindManyP< anab::Calorimetry  > fmcal( recotracks, ev, fRecoTrackCalorimetryTag );
            art::FindManyP< anab::ParticleID   > fmpid( recotracks, ev, fRecoTrackParticleIDTag );
            art::FindManyP< recob::Hit         > fmhit( recotracks, ev, fRecoTrackTag );

            // Loop over tracks associated with primary PFParticles
            for(size_t i = 0; i < trk_assn.size(); ++i) {
              // Get the track-based variables
              std::vector< art::Ptr<anab::Calorimetry> > cal_assn = fmcal.at(trk_assn[i]->ID());
              std::vector< art::Ptr<anab::ParticleID> >  pid_assn = fmpid.at(trk_assn[i]->ID());
              std::vector< art::Ptr<recob::Hit> >        hit_assn = fmhit.at(trk_assn[i]->ID());

              float track_vtx_x = trk_assn[i]->Vertex().X();
              float track_vtx_y = trk_assn[i]->Vertex().Y();
              float track_vtx_z = trk_assn[i]->Vertex().Z();
              float track_end_x = trk_assn[i]->End().X();
              float track_end_y = trk_assn[i]->End().Y();
              float track_end_z = trk_assn[i]->End().Z();

              TVector3 start(track_vtx_x, track_vtx_y, track_vtx_z);
              TVector3 end(track_end_x, track_end_y, track_end_z);
              TVector3 direction = (end - start) * (1./(end - start).Mag());

              // The border for contained tracks should be the edge of the active volume,
              // since this is where we can measure energy up to
              // Find out if one end of a track escapes (if so, MCS)
              bool does_vtx_escape(false);
              for(unsigned int b = 0; b < minx.size(); ++b){
                if(  (track_vtx_x > minx[b])
                    || (track_vtx_x < minx[b])
                    || (track_vtx_y > miny[b])
                    || (track_vtx_y < miny[b])
                    || (track_vtx_z > minz[b])
                    || (track_vtx_z < minz[b])) does_vtx_escape = true;
              }

              bool does_end_escape(false);
              for(unsigned int b = 0; b < minx.size(); ++b){
                if(  (track_end_x > minx[b])
                    || (track_end_x < minx[b])
                    || (track_end_y > miny[b])
                    || (track_end_y < miny[b])
                    || (track_end_z > minz[b])
                    || (track_end_z < minz[b])) does_end_escape = true;
              }

              bool one_end_escapes = true;
              if(does_vtx_escape && does_end_escape)   one_end_escapes = false;
              if(!does_vtx_escape && !does_end_escape) one_end_escapes = false;

              // Loop over PID association
              for ( size_t j = 0; j < pid_assn.size(); ++j ){

                if (!pid_assn[j]) continue;
                if (!pid_assn[j]->PlaneID().isValid) continue;

                // Get the plane number
                int planenum = pid_assn[j]->PlaneID().Plane;

                // Only look at the collection plane, since this is where the dEdx
                // is acquired and we need this for the PIDA values
                if (planenum!=2) continue;

                // Loop over cal association
                for ( size_t k = 0; k < cal_assn.size(); ++k ){

                  if (!cal_assn[k]) continue;
                  if (!cal_assn[k]->PlaneID().isValid) continue;

                  // Get the plane number
                  int planenumcal = cal_assn[k]->PlaneID().Plane;

                  // Only look at the collection plane, since this is where the dEdx
                  // is acquired and we need this for the PIDA values
                  if (planenumcal!=2) continue;

                  // Get associated MCParticle ID using 3 different methods:
                  //    Which particle contributes the most energy to all the hits
                  //    Which particle contributes the reco charge to all the hits
                  //    Which particle is the biggest contributor to all the hitsi
                  fsrp.mc_id_energy = util::TrueParticleIDFromTotalTrueEnergy(hit_assn,fProviderManager);
                  fsrp.mc_id_charge = util::TrueParticleIDFromTotalRecoCharge(hit_assn,fProviderManager);
                  fsrp.mc_id_hits   = util::TrueParticleIDFromTotalRecoHits(hit_assn,fProviderManager);
                  fsrp.hits         = hit_assn.size();
                  fsrp.istrack      = true;
                  fsrp.isshower     = false;
                  
                  /*
                   *      VARIABLES TO FILL
                   */
                  fsrp.chi2_proton         = pid_assn[j]->Chi2Proton();
                  fsrp.chi2_muon           = pid_assn[j]->Chi2Muon();
                  fsrp.chi2_pion           = pid_assn[j]->Chi2Pion();
                  fsrp.pida                = pid_assn[j]->PIDA();
                  fsrp.missing_energy      = pid_assn[j]->MissingE();
                  fsrp.kinetic_energy      = cal_assn[k]->KineticEnergy();
                  fsrp.range               = cal_assn[k]->Range();
                  fsrp.length              = trk_assn[i]->Length();
                  fsrp.vertex              = start;
                  fsrp.end                 = end;
                  fsrp.direction           = direction; 

                  fsrp.dedx_size      = cal_assn[k]->dEdx().size();
                  fsrp.res_range_size = cal_assn[k]->ResidualRange().size();
                  fsrp.pitch_size     = cal_assn[k]->TrkPitchVec().size();
                  for(int l = 0; l < fsrp.dedx_size; ++l)      fsrp.dedx[l]           = cal_assn[k]->dEdx()[l];
                  for(int l = 0; l < fsrp.res_range_size; ++l) fsrp.res_range[l]      = cal_assn[k]->ResidualRange()[l];
                  for(int l = 0; l < fsrp.pitch_size; ++l)     fsrp.pitch[l]          = cal_assn[k]->TrkPitchVec()[l];
/*                  for(int l = fsrp.dedx_size; l < 100000; ++l) fsrp.dedx[l]           = Event::kUnfilled;
                  for(int l = fsrp.res_range_size; l < 100000; ++l) fsrp.res_range[l] = Event::kUnfilled;
                  for(int l = fsrp.pitch_size; l < 100000; ++l)     fsrp.pitch[l]     = Event::kUnfilled;*/

                  // Momentum
                  //    Assign a momentum of 0 to every parameter and then 
                  //    fill with the relevant value based on:
                  //    Whether the track escapes
                  //      - can have an MCS momentum
                  //    If it is contained
                  //      - proton range-based
                  //      - muon range-based
                  //
                  // Further down the line, once PID has been performed properly, 
                  // this should be strictly selected.
                  fsrp.range_momentum_proton = 0.;
                  fsrp.range_momentum_muon   = 0.;
                  fsrp.mcs_momentum_muon     = 0.;
                  double reco_momentum_muon, reco_momentum_proton;
                  if(one_end_escapes){
                    recob::MCSFitResult mcs_result = fMCSFitter->fitMcs(*trk_assn[i]);
                    reco_momentum_muon = mcs_result.bestMomentum();
                    if(reco_momentum_muon < 0) fsrp.mcs_momentum_muon = 0.;
                    fsrp.mcs_momentum_muon = reco_momentum_muon;
                  }
                  else if(!does_vtx_escape && !does_end_escape){
                    reco_momentum_muon   = fRangeFitter.GetTrackMomentum(fsrp.length, 13);
                    reco_momentum_proton = fRangeFitter.GetTrackMomentum(fsrp.length, 2212);
                    if(reco_momentum_muon   == 0) fsrp.range_momentum_muon   = 0.;
                    if(reco_momentum_proton == 0) fsrp.range_momentum_proton = 0.;
                    fsrp.range_momentum_muon   = reco_momentum_muon;
                    fsrp.range_momentum_proton = reco_momentum_proton;
                  }
                } // calo
              } // PID
              reconstructed_interaction.recofinalstate.push_back(fsrp);
            } // tracks
          } // track check

          // Get the shower associations
          //std::vector< art::Ptr<recob::Shower> > shw_assn = fmshw.at(pfparticle.Self());

          // If the showers exist get the information associated with them
          if(showersize) {
            art::FindManyP< recob::Hit > fmhit(recoshowers, ev, fRecoShowerTag);

            // Loop over showers associated with primary PFParticles
            for(int i = 0; i < showersize; ++i){
              // Get the current PFParticle
              auto const shower = recoshowers->at(i);
              std::vector< art::Ptr<recob::Hit> > hit_assn = fmhit.at(shower.ID());

              // Distance s of the shower start from the reconstructed 
              // neutrino vertex
              double s_vtx = sqrt(pow(shower.ShowerStart()[0] - vertex[0],2) + pow(shower.ShowerStart()[1] - vertex[1],2) + pow(shower.ShowerStart()[2] - vertex[2],2));
              int bp = shower.best_plane();

              // Cut of 40 cm for showers to accommodate photon conversion after
              // neutral pion decay
              if(s_vtx < 40) {
                std::vector< art::Ptr<recob::Hit> > hit_assn = fmhit.at(shower.ID());

                // Get associated MCParticle ID using 3 different methods:
                //    Which particle contributes the most energy to all the hits
                //    Which particle contributes the reco charge to all the hits
                //    Which particle is the biggest contributor to all the hits
                fsrp.mc_id_energy = util::TrueParticleIDFromTotalTrueEnergy(hit_assn,fProviderManager);
                fsrp.mc_id_charge = util::TrueParticleIDFromTotalRecoCharge(hit_assn,fProviderManager);
                fsrp.mc_id_hits   = util::TrueParticleIDFromTotalRecoHits(hit_assn,fProviderManager);
                fsrp.hits         = hit_assn.size();
                fsrp.istrack      = false;
                fsrp.isshower     = true;

                fsrp.dedx_size = shower.dEdx().size();
                for(int l = 0; l < fsrp.dedx_size; ++l) fsrp.dedx[l]      = shower.dEdx()[l];
                //for(int l = fsrp.dedx_size; l < 100000; ++l) fsrp.dedx[l] = Event::kUnfilled;
                fsrp.energy        = shower.Energy()[bp];
                fsrp.vertex[0]     = shower.ShowerStart()[0];
                fsrp.vertex[1]     = shower.ShowerStart()[1];
                fsrp.vertex[2]     = shower.ShowerStart()[2];
                fsrp.direction[0]  = shower.Direction()[0];
                fsrp.direction[1]  = shower.Direction()[1];
                fsrp.direction[2]  = shower.Direction()[2];
                fsrp.length        = shower.Length();
                fsrp.opening_angle = shower.OpenAngle();

              } // Vertex check
              reconstructed_interaction.recofinalstate.push_back(fsrp);
            } // showers
          } // shower check
        } // daughters of the neutrino
      } // pfparticles
      fEvent->reco.push_back(reconstructed_interaction);
    } // neutrinos
  } // BuildEventTree
} // namespace core
