#ifndef SBN_ANALYSIS_CALORIMETRYANALYSIS_CXX
#define SBN_ANALYSIS_CALORIMETRYANALYSIS_CXX

#include <iostream>

#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"


#include "lardataobj/AnalysisBase/T0.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbncode/LArRecoProducer/Products/RangeP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"


// #include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "sbncode/CAFMaker/RecoUtils/RecoUtils.h"


namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       CalorimetryAnalysis
// File:        CalorimetryAnalysis.cc
//
//              A basic analysis example
//
// Configuration parameters:
//
// TBD
//
// Created by Giuseppe Cerati (cerati@fnal.gov) on 03/15/2019
//
////////////////////////////////////////////////////////////////////////

// declare helper functions
void TrkDirectionAtXYZ(const recob::Track trk, const double x, const double y, const double z, float out[3]);
void True2RecoMappingXYZ(float t, float x, float y, float z, float out[3]);
float YZtoPlanecoordinate(const float y, const float z, const int plane);
float distance3d(const float& x1, const float& y1, const float& z1,
                  const float& x2, const float& y2, const float& z2);

class CalorimetryAnalysis : public art::EDAnalyzer
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  CalorimetryAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~CalorimetryAnalysis(){};

  /**
     * @brief Analysis function
     */
  void analyze(art::Event const &e) override;

  /**
     * @brief Fill Default info for event associated to neutrino
     */
  void fillDefault();

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree);

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree);


private:

  // function that given a track and its calo info fills NTuple variables
  void FillCalorimetry(art::Event const &e,
            const art::Ptr<recob::PFParticle> &pfp,
            const art::Ptr<recob::Track> &trk,
            const std::vector<art::Ptr<anab::Calorimetry>> &calos,
            const std::vector<art::Ptr<anab::ParticleID>> &pid,
            const std::vector<art::Ptr<recob::Cluster>> &clusters,
            const art::FindManyP<recob::Hit> fmClustertoHit,
            const recob::MCSFitResult &muon_mcs,
            const sbn::RangeP &muon_range,
            const sbn::RangeP &proton_range,
            const bool fData,
            const bool fShrFit,
            const std::vector<std::pair<int, float>> &particleMatches,
            const art::Ptr<simb::MCParticle> &trueParticle,
            const std::vector<const sim::IDE*> &trueParticleIDEs,
            const std::vector<geo::BoxBoundedGeo> &activeVolumes);


  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fPFPproducer;
  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  art::InputTag fT0producer;

  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;
  bool fBacktrack; // do the backtracking needed for this module?

  bool fShrFit; // use shower track-fitter info?

  std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]

  bool fGetCaloID; // get the index of the calorimetry object manually. Needs to be true unless object produced by Pandora hierarchy
  // (This is at least what I foudn empirically)

  TTree* _calo_tree;

  int _run, _sub, _evt;
  // backtracking information
  int _backtracked_pdg;            // PDG code of backtracked particle
  float _backtracked_e;            // energy of backtracked particle
  float _backtracked_purity;       // purity of backtracking
  float _backtracked_completeness; // completeness of backtracking
  float _backtracked_overlay_purity; // purity of overlay

  float _backtracked_px;
  float _backtracked_py;
  float _backtracked_pz;

  float _backtracked_start_x;
  float _backtracked_start_y;
  float _backtracked_start_z;
  float _backtracked_start_t;
  float _backtracked_start_U;
  float _backtracked_start_V;
  float _backtracked_start_Y;
  float _backtracked_sce_start_x;
  float _backtracked_sce_start_y;
  float _backtracked_sce_start_z;
  float _backtracked_sce_start_U;
  float _backtracked_sce_start_V;
  float _backtracked_sce_start_Y;

  bool _backtracked_process_is_stopping;
  bool _backtracked_end_in_tpc;

  uint _generation;    // generation, 1 is primary
  uint _shr_daughters; // number of shower daughters
  uint _trk_daughters; // number of track daughters

  // track information
  int _nplanehits_U;
  int _nplanehits_V;
  int _nplanehits_Y;
  float _trk_score;

  float _trk_theta;
  float _trk_phi;

  float _trk_dir_x;
  float _trk_dir_y;
  float _trk_dir_z;

  float _trk_len;

  float _trk_start_x;
  float _trk_start_y;
  float _trk_start_z;

  float _trk_sce_start_x;
  float _trk_sce_start_y;
  float _trk_sce_start_z;

  float _trk_end_x;
  float _trk_end_y;
  float _trk_end_z;

  float _trk_sce_end_x;
  float _trk_sce_end_y;
  float _trk_sce_end_z;

  float _trk_bragg_p_y;
  float _trk_bragg_mu_y;
  float _trk_bragg_mip_y;
  float _trk_pid_chipr_y;
  float _trk_pid_chika_y;
  float _trk_pid_chipi_y;
  float _trk_pid_chimu_y;
  float _trk_pida_y;

  float _trk_bragg_p_u;
  float _trk_bragg_mu_u;
  float _trk_bragg_mip_u;
  float _trk_pid_chipr_u;
  float _trk_pid_chika_u;
  float _trk_pid_chipi_u;
  float _trk_pid_chimu_u;
  float _trk_pida_u;

  float _trk_bragg_p_v;
  float _trk_bragg_mu_v;
  float _trk_bragg_mip_v;
  float _trk_pid_chipr_v;
  float _trk_pid_chika_v;
  float _trk_pid_chipi_v;
  float _trk_pid_chimu_v;
  float _trk_pida_v;

  float _trk_bragg_p_three_planes;

  float _trk_mcs_muon_mom;
  float _trk_energy_proton;
  float _trk_energy_muon;

  int _longest; // longest track in slice?

  // dedx vector
  std::vector<float> _dqdx_u;
  std::vector<float> _dqdx_v;
  std::vector<float> _dqdx_y;

  std::vector<float> _dedx_u;
  std::vector<float> _dedx_v;
  std::vector<float> _dedx_y;

  std::vector<float> _rr_u;
  std::vector<float> _rr_v;
  std::vector<float> _rr_y;

  std::vector<float> _pitch_u;
  std::vector<float> _pitch_v;
  std::vector<float> _pitch_y;

  std::vector<float> _x_u;
  std::vector<float> _x_v;
  std::vector<float> _x_y;

  std::vector<float> _y_u;
  std::vector<float> _y_v;
  std::vector<float> _y_y;

  std::vector<float> _z_u;
  std::vector<float> _z_v;
  std::vector<float> _z_y;

  std::vector<float> _dir_x_u;
  std::vector<float> _dir_x_v;
  std::vector<float> _dir_x_y;

  std::vector<float> _dir_y_u;
  std::vector<float> _dir_y_v;
  std::vector<float> _dir_y_y;

  std::vector<float> _dir_z_u;
  std::vector<float> _dir_z_v;
  std::vector<float> _dir_z_y;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
CalorimetryAnalysis::CalorimetryAnalysis(const fhicl::ParameterSet &p)
  : art::EDAnalyzer{p}

{
  fPFPproducer  = p.get< art::InputTag > ("PFPproducer","pandora");
  fCALOproducer = p.get< art::InputTag > ("CALOproducer");
  fPIDproducer  = p.get< art::InputTag > ("PIDproducer" );
  fTRKproducer  = p.get< art::InputTag > ("TRKproducer" );
  fT0producer  = p.get< art::InputTag > ("T0producer", "" );

  fBacktrack = p.get<bool>("Backtrack", true);
  fShrFit    = p.get<bool>("ShrFit"   , false);
  fGetCaloID = p.get<bool>("GetCaloID", false);

  fADCtoE = p.get<std::vector<float>>("ADCtoE");

  art::ServiceHandle<art::TFileService> tfs;

  _calo_tree = tfs->make<TTree>("CalorimetryAnalyzer", "Calo Tree");

  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes 
  for (auto const &cryo: geometry->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator iTPC = geometry->begin_TPC(cryo.ID()),
                                    tend = geometry->end_TPC(cryo.ID());
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    while (iTPC != tend) {
      geo::TPCGeo const& TPC = *iTPC;
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
      iTPC++;
    }
     fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: fTPCVolumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }

  setBranches(NULL);
}

void CalorimetryAnalysis::analyze(art::Event const &e)
{
  bool fData = e.isRealData();
  (void) fData;

  _evt = e.event();
  _sub = e.subRun();
  _run = e.run();
  std::cout << "[CalorimetryAnalysis::analyzeEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: "<< _evt << std::endl;

  // Build larpandora info:
  //lar_pandora::LArPandoraHelper larpandora;
  //lar_pandora::PFParticleVector pfparticles;
  //lar_pandora::PFParticleMap particleMap;
  //larpandora.CollectPFParticles(e, "pandora", pfparticles);
  //larpandora.BuildPFParticleMap(pfparticles, particleMap);

  // true particles
  art::ValidHandle<std::vector<simb::MCParticle>> particleHandle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  std::vector<art::Ptr<simb::MCParticle>> particleList;
  art::fill_ptr_vector(particleList, particleHandle);

  art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  art::fill_ptr_vector(PFParticleList, pfparticles);

  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTRKproducer); 

  art::FindManyP<recob::Track> fmTracks(pfparticles, e, fTRKproducer);
  art::FindManyP<recob::Cluster> fmClusters(pfparticles, e, fPFPproducer);
  // art::FindManyP<recob::SpacePoint> fmSpacepoints(pfparticles, e, fPFPproducer); 
  // also get hits...
  art::FindManyP<recob::Hit> fmHits(tracks, e, fTRKproducer);
  art::FindManyP<anab::Calorimetry> fmCalo(tracks, e, fCALOproducer);
  art::FindManyP<anab::ParticleID> fmPID(tracks, e, fPIDproducer);

  art::InputTag mcs_muon  {"pandoraTrackMCS", "muon"};
  art::InputTag range_muon {"pandoraTrackRange", "muon"};
  art::InputTag range_proton {"pandoraTrackRange", "proton"};

  art::FindManyP<recob::MCSFitResult> fmMuonMCS(tracks, e, mcs_muon);
  art::FindManyP<sbn::RangeP> fmMuonRange(tracks, e, range_muon);
  art::FindManyP<sbn::RangeP> fmProtonRange(tracks, e, range_proton);

  for (art::Ptr<recob::PFParticle> p_pfp: PFParticleList) {
    const recob::PFParticle &pfp = *p_pfp;

    fillDefault();

    // grab associated tracks
    const std::vector<art::Ptr<recob::Track>> thisTrack = fmTracks.at(pfp.Self());
    if (thisTrack.size() != 1)
      continue;

    const recob::Track &track = *thisTrack.at(0);

    // get all the data

    // particle data
    const std::vector<art::Ptr<recob::Cluster>> &clusters = fmClusters.at(pfp.Self());

    // track data
    const std::vector<art::Ptr<anab::Calorimetry>> &calo = fmCalo.at(track.ID());
    const std::vector<art::Ptr<anab::ParticleID>> &pid = fmPID.at(track.ID());
    const recob::MCSFitResult &muon_mcs = *fmMuonMCS.at(track.ID()).at(0);
    const sbn::RangeP &muon_range = *fmMuonRange.at(track.ID()).at(0);
    const sbn::RangeP &proton_range = *fmProtonRange.at(track.ID()).at(0);
    
    art::FindManyP<recob::Hit> fmClustertoHit(clusters, e, fPFPproducer);

    // Get the true matching MC particle
    const std::vector<art::Ptr<recob::Hit>> &hits = fmHits.at(track.ID());
    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(hits, true);
    // get the match with the most energy
    int match_id = matches.size() ? std::max_element(matches.begin(), matches.end(),
                                                     [](const auto &a, const auto &b) { return a.second < b.second; })->first : -1;
    float matchE = matches.size() ? std::max_element(matches.begin(), matches.end(),
                                                     [](const auto &a, const auto &b) { return a.second < b.second; })->second : -1;
    art::Ptr<simb::MCParticle> true_particle;
    for (art::Ptr<simb::MCParticle> p: particleList) {
      if (p->TrackId() == match_id) {
        true_particle = p;
        break;
      }
    }
    float trueTrackE = 0.;
    for (auto pair: matches) trueTrackE += pair.second;

    // TODO: fix??
    //
    // ignore tracks with no match
    if (!true_particle) continue;

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    std::vector<const sim::IDE*> particle_ides(bt_serv->TrackIdToSimIDEs_Ps(true_particle->TrackId()));
    FillCalorimetry(e,
          p_pfp,
	  thisTrack.at(0),
          calo,
          pid,
          clusters,
          fmClustertoHit,
          muon_mcs,
          muon_range,
          proton_range,
          false,
          false,
          matches,
          true_particle,
          particle_ides,          
          fActiveVolumes);

  }// for all PFParticles
}// analyzeEvent

void CalorimetryAnalysis::fillDefault()
{
  // backtracking information
  _backtracked_pdg = std::numeric_limits<int>::lowest();            // PDG code of backtracked particle
  _backtracked_e = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_purity = std::numeric_limits<float>::lowest();       // purity of backtracking
  _backtracked_completeness = std::numeric_limits<float>::lowest(); // completeness of backtracking
  _backtracked_overlay_purity = std::numeric_limits<float>::lowest(); // purity of overlay
  _backtracked_process_is_stopping = false;
  _backtracked_end_in_tpc = false;

  _backtracked_px = std::numeric_limits<float>::lowest();
  _backtracked_py = std::numeric_limits<float>::lowest();
  _backtracked_pz = std::numeric_limits<float>::lowest();

  _backtracked_start_x = std::numeric_limits<float>::lowest();
  _backtracked_start_y = std::numeric_limits<float>::lowest();
  _backtracked_start_z = std::numeric_limits<float>::lowest();
  _backtracked_start_t = std::numeric_limits<float>::lowest();
  _backtracked_start_U = std::numeric_limits<float>::lowest();
  _backtracked_start_V = std::numeric_limits<float>::lowest();
  _backtracked_start_Y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_x = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_y = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_z = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_U = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_V = std::numeric_limits<float>::lowest();
  _backtracked_sce_start_Y = std::numeric_limits<float>::lowest();

  _generation = std::numeric_limits<uint>::lowest();
  _shr_daughters = std::numeric_limits<uint>::lowest();
  _trk_daughters = std::numeric_limits<uint>::lowest();

  // track information
  _nplanehits_U = std::numeric_limits<int>::lowest();
  _nplanehits_V = std::numeric_limits<int>::lowest();
  _nplanehits_Y = std::numeric_limits<int>::lowest();
  _trk_score = std::numeric_limits<float>::lowest();

  _trk_theta = std::numeric_limits<float>::lowest();
  _trk_phi = std::numeric_limits<float>::lowest();
  _trk_len = std::numeric_limits<float>::lowest();

  _trk_dir_x = std::numeric_limits<float>::lowest();
  _trk_dir_y = std::numeric_limits<float>::lowest();
  _trk_dir_z = std::numeric_limits<float>::lowest();

  _trk_start_x = std::numeric_limits<float>::lowest();
  _trk_start_y = std::numeric_limits<float>::lowest();
  _trk_start_z = std::numeric_limits<float>::lowest();

  _trk_sce_start_x = std::numeric_limits<float>::lowest();
  _trk_sce_start_y = std::numeric_limits<float>::lowest();
  _trk_sce_start_z = std::numeric_limits<float>::lowest();

  _trk_end_x = std::numeric_limits<float>::lowest();
  _trk_end_y = std::numeric_limits<float>::lowest();
  _trk_end_z = std::numeric_limits<float>::lowest();

  _trk_sce_end_x = std::numeric_limits<float>::lowest();
  _trk_sce_end_y = std::numeric_limits<float>::lowest();
  _trk_sce_end_z = std::numeric_limits<float>::lowest();

  _trk_bragg_p_y = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_y = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_y = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_y = std::numeric_limits<float>::lowest();
  _trk_pid_chika_y = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_y = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_y = std::numeric_limits<float>::lowest();
  _trk_pida_y = std::numeric_limits<float>::lowest();

  _trk_bragg_p_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_u = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_u = std::numeric_limits<float>::lowest();
  _trk_pid_chika_u = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_u = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_u = std::numeric_limits<float>::lowest();
  _trk_pida_u = std::numeric_limits<float>::lowest();

  _trk_bragg_p_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mu_v = std::numeric_limits<float>::lowest();
  _trk_bragg_mip_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipr_v = std::numeric_limits<float>::lowest();
  _trk_pid_chika_v = std::numeric_limits<float>::lowest();
  _trk_pid_chipi_v = std::numeric_limits<float>::lowest();
  _trk_pid_chimu_v = std::numeric_limits<float>::lowest();
  _trk_pida_v = std::numeric_limits<float>::lowest();

  _trk_bragg_p_three_planes = std::numeric_limits<float>::lowest();

  _longest = std::numeric_limits<int>::lowest();

  _trk_mcs_muon_mom = std::numeric_limits<float>::lowest();
  _trk_energy_proton = std::numeric_limits<float>::lowest();
  _trk_energy_muon = std::numeric_limits<float>::lowest();

  // dedx vector
  _dqdx_u.clear();
  _dqdx_v.clear();
  _dqdx_y.clear();

  _dedx_u.clear();
  _dedx_v.clear();
  _dedx_y.clear();

  _rr_u.clear();
  _rr_v.clear();
  _rr_y.clear();

  _pitch_u.clear();
  _pitch_v.clear();
  _pitch_y.clear();

  _x_u.clear();
  _x_v.clear();
  _x_y.clear();

  _y_u.clear();
  _y_v.clear();
  _y_y.clear();

  _z_u.clear();
  _z_v.clear();
  _z_y.clear();

  _dir_x_u.clear();
  _dir_x_v.clear();
  _dir_x_y.clear();

  _dir_y_u.clear();
  _dir_y_v.clear();
  _dir_y_y.clear();

  _dir_z_u.clear();
  _dir_z_v.clear();
  _dir_z_y.clear();
}

void CalorimetryAnalysis::setBranches(TTree *_tree)
{
  _calo_tree->Branch("run", &_run, "run/i");
  _calo_tree->Branch("sub", &_sub, "sub/i");
  _calo_tree->Branch("evt", &_evt, "evt/i");
  // backtracking information
  _calo_tree->Branch("backtracked_pdg", &_backtracked_pdg, "backtracked_pdg/I");            // PDG code of backtracked particle
  _calo_tree->Branch("backtracked_e", &_backtracked_e, "backtracked_e/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_purity", &_backtracked_purity, "backtracked_purity/F");       // purity of backtracking
  _calo_tree->Branch("backtracked_completeness", &_backtracked_completeness, "backtracked_completeness/F"); // completeness of backtracking
  _calo_tree->Branch("backtracked_overlay_purity", &_backtracked_overlay_purity, "backtracked_overlay_purity/F"); // purity of overlay

  _calo_tree->Branch("backtracked_start_x", &_backtracked_start_x, "backtracked_start_x/F");
  _calo_tree->Branch("backtracked_start_y", &_backtracked_start_y, "backtracked_start_y/F");
  _calo_tree->Branch("backtracked_start_z", &_backtracked_start_z, "backtracked_start_z/F");
  _calo_tree->Branch("backtracked_start_t", &_backtracked_start_t, "backtracked_start_t/F");
  // _calo_tree->Branch("backtracked_start_U", &_backtracked_start_U, "backtracked_start_U/F");
  // _calo_tree->Branch("backtracked_start_V", &_backtracked_start_V, "backtracked_start_V/F");
  // _calo_tree->Branch("backtracked_start_Y", &_backtracked_start_Y, "backtracked_start_Y/F");
  _calo_tree->Branch("backtracked_sce_start_x", &_backtracked_sce_start_x, "backtracked_sce_start_x/F");
  _calo_tree->Branch("backtracked_sce_start_y", &_backtracked_sce_start_y, "backtracked_sce_start_y/F");
  _calo_tree->Branch("backtracked_sce_start_z", &_backtracked_sce_start_z, "backtracked_sce_start_z/F");
  // _calo_tree->Branch("backtracked_sce_start_U", &_backtracked_sce_start_U, "backtracked_sce_start_U/F");
  // _calo_tree->Branch("backtracked_sce_start_V", &_backtracked_sce_start_V, "backtracked_sce_start_V/F");
  // _calo_tree->Branch("backtracked_sce_start_Y", &_backtracked_sce_start_Y, "backtracked_sce_start_Y/F");

  _calo_tree->Branch("backtracked_process_is_stopping", &_backtracked_process_is_stopping, "backtracked_process_is_stopping/O");
  _calo_tree->Branch("backtracked_end_in_tpc", &_backtracked_end_in_tpc, "backtracked_end_in_tpc/O");

  _calo_tree->Branch("generation", &_generation, "generation/i");
  _calo_tree->Branch("trk_daughters", &_trk_daughters, "trk_daughters/i");
  _calo_tree->Branch("shr_daughters", &_shr_daughters, "shr_daughters/i");

  // track information
  _calo_tree->Branch("nplanehits_U", &_nplanehits_U, "nplanehits_U/I");
  _calo_tree->Branch("nplanehits_V", &_nplanehits_V, "nplanehits_V/I");
  _calo_tree->Branch("nplanehits_Y", &_nplanehits_Y, "nplanehits_Y/I");
  _calo_tree->Branch("trk_score", &_trk_score, "trk_score/F");

  _calo_tree->Branch("trk_theta", &_trk_theta, "trk_theta/F");
  _calo_tree->Branch("trk_phi", &_trk_phi, "trk_phi/F");
  _calo_tree->Branch("trk_len", &_trk_len, "trk_len/F");

  _calo_tree->Branch("trk_dir_x", &_trk_dir_x, "trk_dir_x/F");
  _calo_tree->Branch("trk_dir_y", &_trk_dir_y, "trk_dir_y/F");
  _calo_tree->Branch("trk_dir_z", &_trk_dir_z, "trk_dir_z/F");

  _calo_tree->Branch("longest", &_longest, "longest/I");

  _calo_tree->Branch("trk_start_x", &_trk_start_x, "trk_start_x/F");
  _calo_tree->Branch("trk_start_y", &_trk_start_y, "trk_start_y/F");
  _calo_tree->Branch("trk_start_z", &_trk_start_z, "trk_start_z/F");

  _calo_tree->Branch("trk_sce_start_x", &_trk_sce_start_x, "trk_sce_start_x/F");
  _calo_tree->Branch("trk_sce_start_y", &_trk_sce_start_y, "trk_sce_start_y/F");
  _calo_tree->Branch("trk_sce_start_z", &_trk_sce_start_z, "trk_sce_start_z/F");

  _calo_tree->Branch("trk_end_x", &_trk_end_x, "trk_end_x/F");
  _calo_tree->Branch("trk_end_y", &_trk_end_y, "trk_end_y/F");
  _calo_tree->Branch("trk_end_z", &_trk_end_z, "trk_end_z/F");

  _calo_tree->Branch("trk_sce_end_x", &_trk_sce_end_x, "trk_sce_end_x/F");
  _calo_tree->Branch("trk_sce_end_y", &_trk_sce_end_y, "trk_sce_end_y/F");
  _calo_tree->Branch("trk_sce_end_z", &_trk_sce_end_z, "trk_sce_end_z/F");

  _calo_tree->Branch("trk_bragg_p_u", &_trk_bragg_p_u, "trk_bragg_p_u/F");
  _calo_tree->Branch("trk_bragg_mu_u", &_trk_bragg_mu_u, "trk_bragg_mu_u/F");
  _calo_tree->Branch("trk_bragg_mip_u", &_trk_bragg_mip_u, "trk_bragg_mip_u/F");
  _calo_tree->Branch("trk_pid_chipr_u", &_trk_pid_chipr_u, "trk_pid_chipr_u/F");
  _calo_tree->Branch("trk_pid_chika_u", &_trk_pid_chika_u, "trk_pid_chika_u/F");
  _calo_tree->Branch("trk_pid_chipi_u", &_trk_pid_chipi_u, "trk_pid_chipi_u/F");
  _calo_tree->Branch("trk_pid_chimu_u", &_trk_pid_chimu_u, "trk_pid_chimu_u/F");
  _calo_tree->Branch("trk_pida_u", &_trk_pida_u, "trk_pida_u/F");

  _calo_tree->Branch("trk_bragg_p_v", &_trk_bragg_p_v, "trk_bragg_p_v/F");
  _calo_tree->Branch("trk_bragg_mu_v", &_trk_bragg_mu_v, "trk_bragg_mu_v/F");
  _calo_tree->Branch("trk_bragg_mip_v", &_trk_bragg_mip_v, "trk_bragg_mip_v/F");
  _calo_tree->Branch("trk_pid_chipr_v", &_trk_pid_chipr_v, "trk_pid_chipr_v/F");
  _calo_tree->Branch("trk_pid_chika_v", &_trk_pid_chika_v, "trk_pid_chika_v/F");
  _calo_tree->Branch("trk_pid_chipi_v", &_trk_pid_chipi_v, "trk_pid_chipi_v/F");
  _calo_tree->Branch("trk_pid_chimu_v", &_trk_pid_chimu_v, "trk_pid_chimu_v/F");
  _calo_tree->Branch("trk_pida_v", &_trk_pida_v, "trk_pida_v/F");

  _calo_tree->Branch("trk_bragg_p_y", &_trk_bragg_p_y, "trk_bragg_p_y/F");
  _calo_tree->Branch("trk_bragg_mu_y", &_trk_bragg_mu_y, "trk_bragg_mu_y/F");
  _calo_tree->Branch("trk_bragg_mip_y", &_trk_bragg_mip_y, "trk_bragg_mip_y/F");
  _calo_tree->Branch("trk_pid_chipr_y", &_trk_pid_chipr_y, "trk_pid_chipr_y/F");
  _calo_tree->Branch("trk_pid_chika_y", &_trk_pid_chika_y, "trk_pid_chika_y/F");
  _calo_tree->Branch("trk_pid_chipi_y", &_trk_pid_chipi_y, "trk_pid_chipi_y/F");
  _calo_tree->Branch("trk_pid_chimu_y", &_trk_pid_chimu_y, "trk_pid_chimu_y/F");
  _calo_tree->Branch("trk_pida_y", &_trk_pida_y, "trk_pida_y/F");

  _calo_tree->Branch("trk_bragg_p_three_planes", &_trk_bragg_p_three_planes, "trk_bragg_p_three_planes/F");

  _calo_tree->Branch("trk_mcs_muon_mom", &_trk_mcs_muon_mom, "trk_mcs_muon_mom/F");
  _calo_tree->Branch("trk_energy_proton", &_trk_energy_proton, "trk_energy_proton/F");
  _calo_tree->Branch("trk_energy_muon", &_trk_energy_muon, "trk_energy_muon/F");

  // dedx vector
  _calo_tree->Branch("dqdx_u", "std::vector<float>", &_dqdx_u);
  _calo_tree->Branch("dqdx_v", "std::vector<float>", &_dqdx_v);
  _calo_tree->Branch("dqdx_y", "std::vector<float>", &_dqdx_y);

  _calo_tree->Branch("dedx_u", "std::vector<float>", &_dedx_u);
  _calo_tree->Branch("dedx_v", "std::vector<float>", &_dedx_v);
  _calo_tree->Branch("dedx_y", "std::vector<float>", &_dedx_y);

  _calo_tree->Branch("rr_u", "std::vector<float>", &_rr_u);
  _calo_tree->Branch("rr_v", "std::vector<float>", &_rr_v);
  _calo_tree->Branch("rr_y", "std::vector<float>", &_rr_y);

  _calo_tree->Branch("pitch_u", "std::vector<float>", &_pitch_u);
  _calo_tree->Branch("pitch_v", "std::vector<float>", &_pitch_v);
  _calo_tree->Branch("pitch_y", "std::vector<float>", &_pitch_y);

  _calo_tree->Branch("x_u", "std::vector<float>", &_x_u);
  _calo_tree->Branch("x_v", "std::vector<float>", &_x_v);
  _calo_tree->Branch("x_y", "std::vector<float>", &_x_y);

  _calo_tree->Branch("y_u", "std::vector<float>", &_y_u);
  _calo_tree->Branch("y_v", "std::vector<float>", &_y_v);
  _calo_tree->Branch("y_y", "std::vector<float>", &_y_y);

  _calo_tree->Branch("z_u", "std::vector<float>", &_z_u);
  _calo_tree->Branch("z_v", "std::vector<float>", &_z_v);
  _calo_tree->Branch("z_y", "std::vector<float>", &_z_y);

  _calo_tree->Branch("dir_x_u", "std::vector<float>", &_dir_x_u);
  _calo_tree->Branch("dir_x_v", "std::vector<float>", &_dir_x_v);
  _calo_tree->Branch("dir_x_y", "std::vector<float>", &_dir_x_y);

  _calo_tree->Branch("dir_y_u", "std::vector<float>", &_dir_y_u);
  _calo_tree->Branch("dir_y_v", "std::vector<float>", &_dir_y_v);
  _calo_tree->Branch("dir_y_y", "std::vector<float>", &_dir_y_y);

  _calo_tree->Branch("dir_z_u", "std::vector<float>", &_dir_z_u);
  _calo_tree->Branch("dir_z_v", "std::vector<float>", &_dir_z_v);
  _calo_tree->Branch("dir_z_y", "std::vector<float>", &_dir_z_y);
}

void CalorimetryAnalysis::resetTTree(TTree *_tree)
{
}


void CalorimetryAnalysis::FillCalorimetry(art::Event const &e,
            const art::Ptr<recob::PFParticle> &pfp,
            const art::Ptr<recob::Track> &trk,
            const std::vector<art::Ptr<anab::Calorimetry>> &calos,
            const std::vector<art::Ptr<anab::ParticleID>> &pids,
            const std::vector<art::Ptr<recob::Cluster>> &clusters,
            const art::FindManyP<recob::Hit> fmClustertoHit,
            const recob::MCSFitResult &muon_mcs,
            const sbn::RangeP &muon_range,
            const sbn::RangeP &proton_range,
            const bool fData,
            const bool fShrFit,
            const std::vector<std::pair<int, float>> &particleMatches,
            const art::Ptr<simb::MCParticle> &trueParticle,
            const std::vector<const sim::IDE*> &trueParticleIDEs,
            const std::vector<geo::BoxBoundedGeo> &activeVolumes)
            //const std::vector<searchingfornues::BtPart> btparts_v,
            //const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart)
{

  // TODO: variables to set
  _generation = -1;
  _shr_daughters = 0;
  _trk_daughters = 0;
  _trk_score = -1;
  _longest = -1;


  std::vector<art::Ptr<recob::Hit>> allHits;

  for (unsigned cluster_i = 0; cluster_i < clusters.size(); cluster_i++) 
  {
    const art::Ptr<recob::Cluster> clus = clusters[cluster_i];
    // get cluster hits
    const std::vector<art::Ptr<recob::Hit>> &cluster_hits = fmClustertoHit.at(cluster_i);
    unsigned nhits = cluster_hits.size();

    _nplanehits_U = 0;
    _nplanehits_V = 0;
    _nplanehits_Y = 0;

    // TODO: does this work in ICARUS???
    if (clus->Plane().Plane == 0)
    {
      _nplanehits_U = nhits;
    }
        else if (clus->Plane().Plane == 1)
    {
      _nplanehits_V = nhits;
    }
        else if (clus->Plane().Plane == 2)
    {
      _nplanehits_Y = nhits;
    }
    for (const auto &hit : cluster_hits)
    {
      allHits.push_back(hit);
    }
  } // for all clusters associated to PFP

  // store Backtracking
  if (!fData)
  {
    _backtracked_e = trueParticle->E();
    _backtracked_pdg = trueParticle->PdgCode();
    float trueTrackE = 0.;
    float matchE = -1;
    for (auto pair: particleMatches) {
      trueTrackE += pair.second;
      if (pair.first == trueParticle->TrackId()) matchE = pair.second;
    }
    assert(matchE > 0.);
    _backtracked_purity = matchE / trueTrackE;

    float trueParticleE = 0.;
    for (auto ide: trueParticleIDEs) trueParticleE += ide->energy;
    _backtracked_completeness = matchE / trueParticleE;

    // TODO: set overlay purity
    

    _backtracked_px = trueParticle->Momentum().Px();
    _backtracked_py = trueParticle->Momentum().Py();
    _backtracked_pz = trueParticle->Momentum().Pz();

    _backtracked_start_x = trueParticle->Position().X();
    _backtracked_start_y = trueParticle->Position().Y();
    _backtracked_start_z = trueParticle->Position().Z();
    _backtracked_start_t = trueParticle->Position().T();

    _backtracked_start_U = YZtoPlanecoordinate(_backtracked_start_y, _backtracked_start_z, 0);
    _backtracked_start_V = YZtoPlanecoordinate(_backtracked_start_y, _backtracked_start_z, 1);
    _backtracked_start_Y = YZtoPlanecoordinate(_backtracked_start_y, _backtracked_start_z, 2);

    // TODO: implement spacecharge
    float reco_st[3] = {_backtracked_start_x, _backtracked_start_y, _backtracked_start_z};
    // TODO: do something special for electrons and photons???
    // if (mcp.pdg == 11 || mcp.pdg == 22)
    // {
    //   reco_st[0] += searchingfornues::x_offset(mcp.start_t);
    // }
    True2RecoMappingXYZ(_backtracked_start_t, _backtracked_start_x, _backtracked_start_y, _backtracked_start_z, reco_st);
    _backtracked_sce_start_x = reco_st[0];
    _backtracked_sce_start_y = reco_st[1];
    _backtracked_sce_start_z = reco_st[2];
    
    _backtracked_sce_start_U = YZtoPlanecoordinate(reco_st[1], reco_st[2], 0);
    _backtracked_sce_start_V = YZtoPlanecoordinate(reco_st[1], reco_st[2], 1);
    _backtracked_sce_start_Y = YZtoPlanecoordinate(reco_st[1], reco_st[2], 2);
    _backtracked_process_is_stopping = (trueParticle->EndProcess() == "Decay" ||
                                        trueParticle->EndProcess() == "CoupledTransportation" ||
                                        trueParticle->EndProcess() == "FastScintillation" ||
                                        trueParticle->EndProcess() == "muMinusCaptureAtRest" ||
                                        trueParticle->EndProcess() == "LArVoxelReadoutScoringProcess");

    _backtracked_end_in_tpc = false;
    for (const geo::BoxBoundedGeo &AV: activeVolumes) {
      if (AV.ContainsPosition(trueParticle->Position().Vect()) && AV.ContainsPosition(trueParticle->EndPosition().Vect())) {
        _backtracked_end_in_tpc = true;
      }
    }
  }

  // Kinetic energy using tabulated stopping power (GeV)
  float mcs_momentum_muon = muon_range.range_p;
  float momentum_proton = proton_range.range_p;
  float energy_proton = std::sqrt(std::pow(momentum_proton, 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
  float energy_muon = std::sqrt(std::pow(mcs_momentum_muon, 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

  _trk_mcs_muon_mom = mcs_momentum_muon;
  _trk_energy_proton = energy_proton;
  _trk_energy_muon = energy_muon;

  _trk_theta = trk->Theta();
  _trk_phi = trk->Phi();
  // TODO: SCE
  _trk_len = trk->Length(); 

  _trk_dir_x = trk->StartDirection().X();
  _trk_dir_y = trk->StartDirection().Y();
  _trk_dir_z = trk->StartDirection().Z();

  _trk_start_x = trk->Start().X();
  _trk_start_y = trk->Start().Y();
  _trk_start_z = trk->Start().Z();

  _trk_end_x = trk->End().X();
  _trk_end_y = trk->End().Y();
  _trk_end_z = trk->End().Z();

  // TODO: SCE
  float _trk_start_sce[3];
  (void) _trk_start_sce;
  _trk_sce_start_x = _trk_start_x;
  _trk_sce_start_y = _trk_start_y;
  _trk_sce_start_z = _trk_start_z;

  float _trk_end_sce[3];
  (void) _trk_end_sce;
  _trk_sce_end_x = _trk_end_x;
  _trk_sce_end_y = _trk_end_y;
  _trk_sce_end_z = _trk_end_z;

  // fill other particle ID
  for (const art::Ptr<anab::ParticleID> pid: pids) {
    auto const &plane = pid->PlaneID().Plane;
    if (plane > 2) continue;

    if (plane == 0) {
       // TODO: bragg variable??
      _trk_bragg_p_u = -1;
      _trk_bragg_mu_u = -1;
      _trk_bragg_mip_u = -1;

      _trk_pid_chipr_u = pid->Chi2Proton();
      _trk_pid_chika_u = pid->Chi2Kaon();
      _trk_pid_chipi_u = pid->Chi2Pion();
      _trk_pid_chimu_u = pid->Chi2Muon();
      _trk_pida_u = pid->PIDA();
    }
    else if (plane == 1) {
       // TODO: bragg variable??
      _trk_bragg_p_v = -1;
      _trk_bragg_mu_v = -1;
      _trk_bragg_mip_v = -1;

      _trk_pid_chipr_v = pid->Chi2Proton();
      _trk_pid_chika_v = pid->Chi2Kaon();
      _trk_pid_chipi_v = pid->Chi2Pion();
      _trk_pid_chimu_v = pid->Chi2Muon();
      _trk_pida_v = pid->PIDA();
    }
    else if (plane == 2) {
       // TODO: bragg variable??
      _trk_bragg_p_y = -1;
      _trk_bragg_mu_y = -1;
      _trk_bragg_mip_y = -1;

      _trk_pid_chipr_y = pid->Chi2Proton();
      _trk_pid_chika_y = pid->Chi2Kaon();
      _trk_pid_chipi_y = pid->Chi2Pion();
      _trk_pid_chimu_y = pid->Chi2Muon();
    }
  }

  // fill Calorimetry
  for (const art::Ptr<anab::Calorimetry> calo : calos)
  {
    // TODO: will this work for ICARUS?
    auto const& plane = calo->PlaneID().Plane;
    if (plane>2)
    {
      continue;
    }
    auto const& xyz_v = calo->XYZ();

    if (plane == 0)
    {
      _dqdx_u = calo->dQdx();
      _rr_u = calo->ResidualRange();
      _pitch_u = calo->TrkPitchVec();
      for (auto xyz : xyz_v)
      {
        _x_u.push_back(xyz.X());
        _y_u.push_back(xyz.Y());
        _z_u.push_back(xyz.Z());

        float _dir_u[3];
        TrkDirectionAtXYZ(*trk, xyz.X(), xyz.Y(), xyz.Z(), _dir_u);
        _dir_x_u.push_back(_dir_u[0]);
        _dir_y_u.push_back(_dir_u[1]);
        _dir_z_u.push_back(_dir_u[2]);
      }
      // TODO: ShrFit???
      // if (fShrFit) { _dedx_u = searchingfornues::GetdEdxfromdQdx(_dqdx_u, _x_u, _y_u, _z_u, 2.1, fADCtoE[plane]); }
      _dedx_u = calo->dEdx();
    }
    else if (plane == 1)
    {
      _dqdx_v = calo->dQdx();
      _rr_v = calo->ResidualRange();
      _pitch_v = calo->TrkPitchVec();
      for (auto xyz : xyz_v)
      {
        _x_v.push_back(xyz.X());
        _y_v.push_back(xyz.Y());
        _z_v.push_back(xyz.Z());

        float _dir_v[3];
        TrkDirectionAtXYZ(*trk, xyz.X(), xyz.Y(), xyz.Z(), _dir_v);
        _dir_x_v.push_back(_dir_v[0]);
        _dir_y_v.push_back(_dir_v[1]);
        _dir_z_v.push_back(_dir_v[2]);
      }
      // TODO: ShrFit???
      // if (fShrFit) { _dedx_v = searchingfornues::GetdEdxfromdQdx(_dqdx_v, _x_v, _y_v, _z_v, 2.1, fADCtoE[plane]); }
      _dedx_v = calo->dEdx();
    }
    else if (plane == 2) //collection
    {
      _dqdx_y = calo->dQdx();
      _rr_y = calo->ResidualRange();
      _pitch_y = calo->TrkPitchVec();
      for (auto xyz : xyz_v)
      {
        _x_y.push_back(xyz.X());
        _y_y.push_back(xyz.Y());
        _z_y.push_back(xyz.Z());

        float _dir_y[3];
        TrkDirectionAtXYZ(*trk, xyz.X(), xyz.Y(), xyz.Z(), _dir_y);
        _dir_x_y.push_back(_dir_y[0]);
        _dir_y_y.push_back(_dir_y[1]);
        _dir_z_y.push_back(_dir_y[2]);
      }
      // TODO: ShrFit???
      // if (fShrFit) { _dedx_y = searchingfornues::GetdEdxfromdQdx(_dqdx_y, _x_y, _y_y, _z_y, 2.1, fADCtoE[plane]); }
      _dedx_y = calo->dEdx();
    }
  }

  _calo_tree->Fill();
}

// Define helper functions
  void TrkDirectionAtXYZ(const recob::Track trk, const double x, const double y, const double z, float out[3])
  {
    float min_dist = 100;
    size_t i_min = -1;
    for(size_t i=0; i < trk.NumberTrajectoryPoints(); i++)
    {
      if (trk.HasValidPoint(i))
      { // check this point is valid
        auto point_i = trk.LocationAtPoint(i);
        float distance = distance3d((double)point_i.X(), (double)point_i.Y(), (double)point_i.Z(),
                  x, y, z);
        if (distance < min_dist)
        {
          min_dist = distance;
          i_min = i;
        }
      }// if point is valid
    }// for all track points

    auto direction = trk.DirectionAtPoint(i_min);
    out[0] = (float)direction.X();
    out[1] = (float)direction.Y();
    out[2] = (float)direction.Z();

    float norm;
    norm = out[0]*out[0] + out[1]*out[1] + out[2]*out[2];
    if (fabs(norm -1) > 0.001)
    {
      std::cout << "i_min = " << i_min << std::endl;
      std::cout << "minimum distance = " << min_dist << std::endl;
      std::cout << "out[0], out[1], out[2] = " << out[0] << " , " << out[1] << " , " << out[2] << std::endl;
      std::cout << "norm = " << norm << std::endl;
    }
  }
// TODO: include space charge
void True2RecoMappingXYZ(float t,float x, float y, float z, float out[3]) {
  (void) t;
  out[0] = x;
  out[1] = y;
  out[2] = z;
}

  float YZtoPlanecoordinate(const float y, const float z, const int plane)
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    double _wire2cm = geom->WirePitch(0, 0, 0);
    return geom->WireCoordinate(y, z, geo::PlaneID(0, 0, plane)) * _wire2cm;
  }

  float distance3d(const float& x1, const float& y1, const float& z1,
                  const float& x2, const float& y2, const float& z2)
  {
    return sqrt((x1-x2)*(x1-x2) +
                (y1-y2)*(y1-y2) +
                (z1-z2)*(z1-z2));

  }
} // namespace analysis

DEFINE_ART_MODULE(analysis::CalorimetryAnalysis)
#endif
