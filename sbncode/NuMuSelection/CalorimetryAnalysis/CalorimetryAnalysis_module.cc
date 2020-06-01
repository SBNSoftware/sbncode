#ifndef SBN_ANALYSIS_CALORIMETRYANALYSIS_CXX
#define SBN_ANALYSIS_CALORIMETRYANALYSIS_CXX

#include <iostream>
#include <map>

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
#include "lardataobj/RecoBase/Wire.h"
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

#include "larcorealg/GeoAlgo/GeoAlgo.h"



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
float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                      const std::vector<geoalgo::AABox> &boxes);
void FillHits(const std::vector<art::Ptr<recob::Hit>> &hits,
              std::vector<float> &charge_u,
              std::vector<float> &charge_v,
              std::vector<float> &charge_y,
              std::vector<unsigned> &channel_u,
              std::vector<unsigned> &channel_v,
              std::vector<unsigned> &channel_y,
              std::vector<int> &multiplicity_u,
              std::vector<int> &multiplicity_v,
              std::vector<int> &multiplicity_y,
              std::vector<float> &width_u,
              std::vector<float> &width_v,
              std::vector<float> &width_y,
              std::vector<float> &time_u,
              std::vector<float> &time_v,
              std::vector<float> &time_y,
              bool use_integral=true);

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

  void respondToOpenInputFile(const art::FileBlock& fb) {
    (void) fb;
    _fileno ++;
  }

private:

  // function that given a track and its calo info fills NTuple variables
  void FillCalorimetry(art::Event const &e,
            const std::vector<art::Ptr<recob::Wire>> &allWires,
            const std::vector<art::Ptr<recob::Hit>> &allHits,
            const art::Ptr<recob::PFParticle> &pfp,
            const art::Ptr<recob::Track> &trk,
            const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
            const std::vector<art::Ptr<anab::Calorimetry>> &calos,
            const std::vector<art::Ptr<anab::ParticleID>> &pid,
            const std::vector<art::Ptr<recob::Hit>> &areaHits,
            const std::vector<art::Ptr<recob::Hit>> &trkHits,
            const std::vector<art::Ptr<recob::Hit>> &caloHits,
            const recob::MCSFitResult &muon_mcs,
            const sbn::RangeP &muon_range,
            const sbn::RangeP &proton_range,
            const bool fData,
            const bool fShrFit,
            const std::vector<std::pair<int, float>> &particleMatches,
            const art::Ptr<simb::MCParticle> &trueParticle,
            const std::array<std::map<unsigned, std::vector<const sim::IDE*>>, 3> particle_ide_map,
            const std::array<std::vector<const sim::IDE*>,3> &allIDEs,
            const std::vector<geo::BoxBoundedGeo> &activeVolumes);


  TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
  TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

  art::InputTag fPFPproducer;
  art::InputTag fCALOproducer;
  art::InputTag fPIDproducer;
  art::InputTag fTRKproducer;
  art::InputTag fT0producer;
  art::InputTag fHitFilterproducer;
  art::InputTag fAreaHitproducer;
  bool fAllTrueEnergyDeposits;

  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;
  bool fBacktrack; // do the backtracking needed for this module?

  bool fShrFit; // use shower track-fitter info?

  std::vector<float> fADCtoE; // vector of ADC to # of e- conversion [to be taken from production reco2 fhicl files]

  bool fGetCaloID; // get the index of the calorimetry object manually. Needs to be true unless object produced by Pandora hierarchy
  // (This is at least what I foudn empirically)

  TTree* _calo_tree;

  int _fileno;
  int _run, _sub, _evt;
  // backtracking information
  int _backtracked_pdg;            // PDG code of backtracked particle
  float _backtracked_e;            // energy of backtracked particle
  float _backtracked_deposited_e_u;            // energy of backtracked particle
  float _backtracked_deposited_e_v;            // energy of backtracked particle
  float _backtracked_deposited_e_y;            // energy of backtracked particle
  float _backtracked_purity;       // purity of backtracking
  float _backtracked_completeness; // completeness of backtracking
  float _backtracked_overlay_purity; // purity of overlay
  float _backtracked_q_u;
  float _backtracked_q_v;
  float _backtracked_q_y;
  float _backtracked_qtotal_u;
  float _backtracked_qtotal_v;
  float _backtracked_qtotal_y;

  std::vector<float> _backtracked_pc_q_u;
  std::vector<float> _backtracked_pc_q_v;
  std::vector<float> _backtracked_pc_q_y;

  std::vector<unsigned> _backtracked_c_u;
  std::vector<unsigned> _backtracked_c_v;
  std::vector<unsigned> _backtracked_c_y;

  float _backtracked_theta;
  float _backtracked_phi;

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

  float _backtracked_length;
  float _backtracked_length_start_to_end;

  bool _backtracked_process_is_stopping;
  bool _backtracked_end_in_tpc;

  uint _generation;    // generation, 1 is primary
  uint _shr_daughters; // number of shower daughters
  uint _trk_daughters; // number of track daughters
  uint _daughters;
  uint _n_pfp;

  // track information
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

  std::vector<float> _integrated_charge_u;
  std::vector<float> _integrated_charge_v;
  std::vector<float> _integrated_charge_y;

  std::vector<float> _allhit_charge_u;
  std::vector<float> _allhit_charge_v;
  std::vector<float> _allhit_charge_y;
  std::vector<unsigned> _allhit_channel_u;
  std::vector<unsigned> _allhit_channel_v;
  std::vector<unsigned> _allhit_channel_y;
  std::vector<int> _allhit_multiplicity_u;
  std::vector<int> _allhit_multiplicity_v;
  std::vector<int> _allhit_multiplicity_y;
  std::vector<float> _allhit_width_u;
  std::vector<float> _allhit_width_v;
  std::vector<float> _allhit_width_y;
  std::vector<float> _allhit_time_u;
  std::vector<float> _allhit_time_v;
  std::vector<float> _allhit_time_y;

  std::vector<float> _trkhit_charge_u;
  std::vector<float> _trkhit_charge_v;
  std::vector<float> _trkhit_charge_y;
  std::vector<unsigned> _trkhit_channel_u;
  std::vector<unsigned> _trkhit_channel_v;
  std::vector<unsigned> _trkhit_channel_y;
  std::vector<int> _trkhit_multiplicity_u;
  std::vector<int> _trkhit_multiplicity_v;
  std::vector<int> _trkhit_multiplicity_y;
  std::vector<float> _trkhit_width_u;
  std::vector<float> _trkhit_width_v;
  std::vector<float> _trkhit_width_y;
  std::vector<float> _trkhit_time_u;
  std::vector<float> _trkhit_time_v;
  std::vector<float> _trkhit_time_y;

  std::vector<float> _calohit_charge_u;
  std::vector<float> _calohit_charge_v;
  std::vector<float> _calohit_charge_y;
  std::vector<unsigned> _calohit_channel_u;
  std::vector<unsigned> _calohit_channel_v;
  std::vector<unsigned> _calohit_channel_y;
  std::vector<int> _calohit_multiplicity_u;
  std::vector<int> _calohit_multiplicity_v;
  std::vector<int> _calohit_multiplicity_y;
  std::vector<float> _calohit_width_u;
  std::vector<float> _calohit_width_v;
  std::vector<float> _calohit_width_y;
  std::vector<float> _calohit_time_u;
  std::vector<float> _calohit_time_v;
  std::vector<float> _calohit_time_y;

  std::vector<float> _sumhit_charge_u;
  std::vector<float> _sumhit_charge_v;
  std::vector<float> _sumhit_charge_y;
  std::vector<unsigned> _sumhit_channel_u;
  std::vector<unsigned> _sumhit_channel_v;
  std::vector<unsigned> _sumhit_channel_y;
  std::vector<int> _sumhit_multiplicity_u;
  std::vector<int> _sumhit_multiplicity_v;
  std::vector<int> _sumhit_multiplicity_y;
  std::vector<float> _sumhit_width_u;
  std::vector<float> _sumhit_width_v;
  std::vector<float> _sumhit_width_y;
  std::vector<float> _sumhit_time_u;
  std::vector<float> _sumhit_time_v;
  std::vector<float> _sumhit_time_y;

  std::vector<float> _areahit_charge_u;
  std::vector<float> _areahit_charge_v;
  std::vector<float> _areahit_charge_y;
  std::vector<unsigned> _areahit_channel_u;
  std::vector<unsigned> _areahit_channel_v;
  std::vector<unsigned> _areahit_channel_y;
  std::vector<int> _areahit_multiplicity_u;
  std::vector<int> _areahit_multiplicity_v;
  std::vector<int> _areahit_multiplicity_y;
  std::vector<float> _areahit_width_u;
  std::vector<float> _areahit_width_v;
  std::vector<float> _areahit_width_y;
  std::vector<float> _areahit_time_u;
  std::vector<float> _areahit_time_v;
  std::vector<float> _areahit_time_y;

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

  std::vector<float> _sx;
  std::vector<float> _sy;
  std::vector<float> _sz;

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
  fAreaHitproducer = p.get<art::InputTag> ("AreaHitproducer", "areahit");
  fAllTrueEnergyDeposits = p.get<bool>("AllTrueEnergyDeposits", true);
  fHitFilterproducer = p.get<art::InputTag>("HitFilterproducer", "filtgoodhit"); 

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

  _fileno = 0;

  setBranches(NULL);
}

void CalorimetryAnalysis::analyze(art::Event const &e)
{
  bool fData = e.isRealData();
  (void) fData;

  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

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

  // recob Hits
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>("gaushit");
  std::vector<art::Ptr<recob::Hit>> hitList;
  art::fill_ptr_vector(hitList, hitHandle);
  // recob Wires
  art::Handle<std::vector<recob::Wire>> wireHandle;
  e.getByLabel("caldata", wireHandle);

  std::vector<art::Ptr<recob::Wire>> wireList;
  if (wireHandle.isValid()) {
    art::fill_ptr_vector(wireList, wireHandle);
  }

  // true particles
  art::ValidHandle<std::vector<simb::MCParticle>> particleHandle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  std::vector<art::Ptr<simb::MCParticle>> particleList;
  art::fill_ptr_vector(particleList, particleHandle);

  art::ValidHandle<std::vector<sim::SimChannel>> simchannelHandle = e.getValidHandle<std::vector<sim::SimChannel>>("largeant");
  std::vector<art::Ptr<sim::SimChannel>> simchannelList;
  art::fill_ptr_vector(simchannelList, simchannelHandle);

  art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  art::fill_ptr_vector(PFParticleList, pfparticles);
  _n_pfp = pfparticles->size();

  art::ValidHandle<std::vector<recob::Track>> tracks = e.getValidHandle<std::vector<recob::Track>>(fTRKproducer); 

  art::FindManyP<recob::SpacePoint> fmSpacePoint(PFParticleList, e, fPFPproducer);

  art::FindManyP<recob::Track> fmTracks(pfparticles, e, fTRKproducer);
  // art::FindManyP<recob::SpacePoint> fmSpacepoints(pfparticles, e, fPFPproducer); 
  // also get hits...
  art::FindManyP<recob::Hit> fmtrkHits(tracks, e, fTRKproducer);
  art::FindManyP<recob::Hit> fmareaHits(tracks, e, fAreaHitproducer);
  art::FindManyP<recob::Hit> fmcaloHits(tracks, e, fHitFilterproducer);

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

    unsigned parent_ind = pfp.Parent();
    std::cout << "Parent ind: " << parent_ind << std::endl;
    if (parent_ind == recob::PFParticle::kPFParticlePrimary || parent_ind >= PFParticleList.size()) {
      std::cout << "Parent doesn't exist.\n";
      continue;
    }
    // check if parent is the primary
    const recob::PFParticle &parent = *PFParticleList.at(pfp.Parent());
    if (!parent.IsPrimary()) {
      std::cout << "Parent not primary\n";
      continue;
    }

    const recob::Track &track = *thisTrack.at(0);

    // get all the data

    // track data
    const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints = fmSpacePoint.at(pfp.Self());
    const std::vector<art::Ptr<anab::Calorimetry>> &calo = fmCalo.at(track.ID());
    const std::vector<art::Ptr<anab::ParticleID>> &pid = fmPID.at(track.ID());
    const recob::MCSFitResult &muon_mcs = *fmMuonMCS.at(track.ID()).at(0);
    const sbn::RangeP &muon_range = *fmMuonRange.at(track.ID()).at(0);
    const sbn::RangeP &proton_range = *fmProtonRange.at(track.ID()).at(0);

    std::vector<art::Ptr<recob::Hit>> emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &trkHits  = fmtrkHits.at(track.ID());
    const std::vector<art::Ptr<recob::Hit>> &areaHits = fmareaHits.isValid() ?  fmareaHits.at(track.ID()) : emptyHitVector;
    const std::vector<art::Ptr<recob::Hit>> &caloHits = fmcaloHits.at(track.ID());
    
    // Get the true matching MC particle
    std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(trkHits, true);
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
    // only take particle that matches to true instigator
    if (true_particle->Process() != "primary") continue;


    std::cout << "Track: " << track.ID() << " len: " << track.Length() << " matched to particle: " << true_particle->TrackId() << " frac: " << (matchE/trueTrackE) << std::endl;

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;

    std::array<std::map<unsigned, std::vector<const sim::IDE*>>, 3> particle_ide_map;
    for (const art::Ptr<sim::SimChannel> sc: simchannelList) {
      for (const auto &pair: sc->TDCIDEMap()) {
        for (const sim::IDE &ide: pair.second) {
          if (ide.trackID == true_particle->TrackId() || fAllTrueEnergyDeposits) {
            unsigned plane_id = geometry->ChannelToWire(sc->Channel()).at(0).Plane;
            particle_ide_map[plane_id][sc->Channel()].push_back(&ide);
          }
        }
      }
    } 

    std::array<std::vector<const sim::IDE*>, 3> allIDEs;
    for (unsigned i = 0; i < particleList.size(); i++) {
      int G4ID = particleList[i]->TrackId();
      for (unsigned plane = 0; plane < 3; plane++) {
        // build the plane ID
        geo::PlaneID plane_id(0, 0, plane);
        std::vector<const sim::IDE*> this_ides = bt_serv->TrackIdToSimIDEs_Ps(G4ID, geometry->View(plane_id));
        allIDEs[plane].insert(allIDEs[plane].end(), this_ides.begin(), this_ides.end());
      }
    }

    FillCalorimetry(e,
          wireList,
          hitList,
          p_pfp,
	  thisTrack.at(0),
          spacepoints,
          calo,
          pid,
          areaHits,
          trkHits,
          caloHits,
          muon_mcs,
          muon_range,
          proton_range,
          false,
          false,
          matches,
          true_particle,
          particle_ide_map, 
          allIDEs,
          fActiveVolumes);

  }// for all PFParticles
}// analyzeEvent

void CalorimetryAnalysis::fillDefault()
{
  // backtracking information
  _backtracked_pdg = std::numeric_limits<int>::lowest();            // PDG code of backtracked particle
  _backtracked_e = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_deposited_e_u = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_deposited_e_v = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_deposited_e_y = std::numeric_limits<float>::lowest();            // energy of backtracked particle
  _backtracked_q_u = std::numeric_limits<float>::lowest();
  _backtracked_q_v = std::numeric_limits<float>::lowest();
  _backtracked_q_y = std::numeric_limits<float>::lowest();
  _backtracked_qtotal_u = std::numeric_limits<float>::lowest();
  _backtracked_qtotal_v = std::numeric_limits<float>::lowest();
  _backtracked_qtotal_y = std::numeric_limits<float>::lowest();
  _backtracked_purity = std::numeric_limits<float>::lowest();       // purity of backtracking
  _backtracked_completeness = std::numeric_limits<float>::lowest(); // completeness of backtracking
  _backtracked_overlay_purity = std::numeric_limits<float>::lowest(); // purity of overlay
  _backtracked_process_is_stopping = false;
  _backtracked_end_in_tpc = false;

  _backtracked_pc_q_u.clear();
  _backtracked_pc_q_v.clear();
  _backtracked_pc_q_y.clear();

  _backtracked_c_u.clear();
  _backtracked_c_v.clear();
  _backtracked_c_y.clear();

  _backtracked_length = std::numeric_limits<float>::lowest();
  _backtracked_length_start_to_end = std::numeric_limits<float>::lowest();

  _backtracked_theta = std::numeric_limits<float>::lowest();
  _backtracked_phi = std::numeric_limits<float>::lowest();

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
  _daughters = std::numeric_limits<uint>::lowest();

  // track information
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


  _integrated_charge_u.clear();
  _integrated_charge_v.clear();
  _integrated_charge_y.clear();

  _trkhit_charge_u.clear();
  _trkhit_charge_v.clear();
  _trkhit_charge_y.clear();
  _trkhit_channel_u.clear();
  _trkhit_channel_v.clear();
  _trkhit_channel_y.clear();
  _trkhit_multiplicity_u.clear();
  _trkhit_multiplicity_v.clear();
  _trkhit_multiplicity_y.clear();
  _trkhit_width_u.clear();
  _trkhit_width_v.clear();
  _trkhit_width_y.clear();
  _trkhit_time_u.clear();
  _trkhit_time_v.clear();
  _trkhit_time_y.clear();

  _areahit_charge_u.clear();
  _areahit_charge_v.clear();
  _areahit_charge_y.clear();
  _areahit_channel_u.clear();
  _areahit_channel_v.clear();
  _areahit_channel_y.clear();
  _areahit_multiplicity_u.clear();
  _areahit_multiplicity_v.clear();
  _areahit_multiplicity_y.clear();
  _areahit_width_u.clear();
  _areahit_width_v.clear();
  _areahit_width_y.clear();
  _areahit_time_u.clear();
  _areahit_time_v.clear();
  _areahit_time_y.clear();

  _allhit_charge_u.clear();
  _allhit_charge_v.clear();
  _allhit_charge_y.clear();
  _allhit_channel_u.clear();
  _allhit_channel_v.clear();
  _allhit_channel_y.clear();
  _allhit_multiplicity_u.clear();
  _allhit_multiplicity_v.clear();
  _allhit_multiplicity_y.clear();
  _allhit_width_u.clear();
  _allhit_width_v.clear();
  _allhit_width_y.clear();
  _allhit_time_u.clear();
  _allhit_time_v.clear();
  _allhit_time_y.clear();

  _calohit_charge_u.clear();
  _calohit_charge_v.clear();
  _calohit_charge_y.clear();
  _calohit_channel_u.clear();
  _calohit_channel_v.clear();
  _calohit_channel_y.clear();
  _calohit_multiplicity_u.clear();
  _calohit_multiplicity_v.clear();
  _calohit_multiplicity_y.clear();
  _calohit_width_u.clear();
  _calohit_width_v.clear();
  _calohit_width_y.clear();
  _calohit_time_u.clear();
  _calohit_time_v.clear();
  _calohit_time_y.clear();

  _sumhit_charge_u.clear();
  _sumhit_charge_v.clear();
  _sumhit_charge_y.clear();
  _sumhit_channel_u.clear();
  _sumhit_channel_v.clear();
  _sumhit_channel_y.clear();
  _sumhit_multiplicity_u.clear();
  _sumhit_multiplicity_v.clear();
  _sumhit_multiplicity_y.clear();
  _sumhit_width_u.clear();
  _sumhit_width_v.clear();
  _sumhit_width_y.clear();
  _sumhit_time_u.clear();
  _sumhit_time_v.clear();
  _sumhit_time_y.clear();

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

  _sx.clear();
  _sy.clear();
  _sz.clear();

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
  _calo_tree->Branch("fileno", &_fileno, "fileno/i");
  _calo_tree->Branch("run", &_run, "run/i");
  _calo_tree->Branch("sub", &_sub, "sub/i");
  _calo_tree->Branch("evt", &_evt, "evt/i");
  // backtracking information
  _calo_tree->Branch("backtracked_pdg", &_backtracked_pdg, "backtracked_pdg/I");            // PDG code of backtracked particle
  _calo_tree->Branch("backtracked_e", &_backtracked_e, "backtracked_e/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_deposited_e_u", &_backtracked_deposited_e_u, "backtracked_deposited_e_u/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_deposited_e_v", &_backtracked_deposited_e_v, "backtracked_deposited_e_v/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_deposited_e_y", &_backtracked_deposited_e_y, "backtracked_deposited_e_y/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_q_u", &_backtracked_q_u, "backtracked_q_u/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_q_v", &_backtracked_q_v, "backtracked_q_v/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_q_y", &_backtracked_q_y, "backtracked_q_y/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_qtotal_u", &_backtracked_qtotal_u, "backtracked_qtotal_u/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_qtotal_v", &_backtracked_qtotal_v, "backtracked_qtotal_v/F");            // energy of backtracked particle
  _calo_tree->Branch("backtracked_qtotal_y", &_backtracked_qtotal_y, "backtracked_qtotal_y/F");            // energy of backtracked particle

  _calo_tree->Branch("backtracked_pc_q_u", "std::vector<float>", &_backtracked_pc_q_u);
  _calo_tree->Branch("backtracked_pc_q_v", "std::vector<float>", &_backtracked_pc_q_v);
  _calo_tree->Branch("backtracked_pc_q_y", "std::vector<float>", &_backtracked_pc_q_y);

  _calo_tree->Branch("backtracked_c_u", "std::vector<unsigned>", &_backtracked_c_u);
  _calo_tree->Branch("backtracked_c_v", "std::vector<unsigned>", &_backtracked_c_v);
  _calo_tree->Branch("backtracked_c_y", "std::vector<unsigned>", &_backtracked_c_y);

  _calo_tree->Branch("backtracked_theta", &_backtracked_theta, "backtracked_theta/F");
  _calo_tree->Branch("backtracked_phi", &_backtracked_phi, "backtracked_phi/F");

  _calo_tree->Branch("backtracked_px", &_backtracked_px, "backtracked_px/F");
  _calo_tree->Branch("backtracked_py", &_backtracked_py, "backtracked_py/F");
  _calo_tree->Branch("backtracked_pz", &_backtracked_pz, "backtracked_pz/F");

  _calo_tree->Branch("backtracked_length", &_backtracked_length,"backtracked_length/F");
  _calo_tree->Branch("backtracked_length_start_to_end", &_backtracked_length_start_to_end, "backtracked_length_start_to_end/F");

  _calo_tree->Branch("backtracked_purity", &_backtracked_purity, "backtracked_purity/F");       // purity of backtracking
  _calo_tree->Branch("backtracked_completeness", &_backtracked_completeness, "backtracked_completeness/F"); // completeness of backtracking
  _calo_tree->Branch("backtracked_overlay_purity", &_backtracked_overlay_purity, "backtracked_overlay_purity/F"); // purity of overlay

  _calo_tree->Branch("integrated_charge_u", "std::vector<float>", &_integrated_charge_u);
  _calo_tree->Branch("integrated_charge_v", "std::vector<float>", &_integrated_charge_v);
  _calo_tree->Branch("integrated_charge_y", "std::vector<float>", &_integrated_charge_y);

  _calo_tree->Branch("areahit_charge_u", "std::vector<float>", &_areahit_charge_u);
  _calo_tree->Branch("areahit_charge_v", "std::vector<float>", &_areahit_charge_v);
  _calo_tree->Branch("areahit_charge_y", "std::vector<float>", &_areahit_charge_y);
  _calo_tree->Branch("areahit_channel_u", "std::vector<unsigned>", &_areahit_channel_u);
  _calo_tree->Branch("areahit_channel_v", "std::vector<unsigned>", &_areahit_channel_v);
  _calo_tree->Branch("areahit_channel_y", "std::vector<unsigned>", &_areahit_channel_y);
  _calo_tree->Branch("areahit_multiplicity_u", "std::vector<int>", &_areahit_multiplicity_u); 
  _calo_tree->Branch("areahit_multiplicity_v", "std::vector<int>", &_areahit_multiplicity_v); 
  _calo_tree->Branch("areahit_multiplicity_y", "std::vector<int>", &_areahit_multiplicity_y); 
  _calo_tree->Branch("areahit_width_u", "std::vector<float>", &_areahit_width_u);
  _calo_tree->Branch("areahit_width_v", "std::vector<float>", &_areahit_width_v);
  _calo_tree->Branch("areahit_width_y", "std::vector<float>", &_areahit_width_y);
  _calo_tree->Branch("areahit_time_u", "std::vector<float>", &_areahit_time_u);
  _calo_tree->Branch("areahit_time_v", "std::vector<float>", &_areahit_time_v);
  _calo_tree->Branch("areahit_time_y", "std::vector<float>", &_areahit_time_y);

  _calo_tree->Branch("allhit_charge_u", "std::vector<float>", &_allhit_charge_u);
  _calo_tree->Branch("allhit_charge_v", "std::vector<float>", &_allhit_charge_v);
  _calo_tree->Branch("allhit_charge_y", "std::vector<float>", &_allhit_charge_y);
  _calo_tree->Branch("allhit_channel_u", "std::vector<unsigned>", &_allhit_channel_u);
  _calo_tree->Branch("allhit_channel_v", "std::vector<unsigned>", &_allhit_channel_v);
  _calo_tree->Branch("allhit_channel_y", "std::vector<unsigned>", &_allhit_channel_y);
  _calo_tree->Branch("allhit_multiplicity_u", "std::vector<int>", &_allhit_multiplicity_u); 
  _calo_tree->Branch("allhit_multiplicity_v", "std::vector<int>", &_allhit_multiplicity_v); 
  _calo_tree->Branch("allhit_multiplicity_y", "std::vector<int>", &_allhit_multiplicity_y); 
  _calo_tree->Branch("allhit_width_u", "std::vector<float>", &_allhit_width_u);
  _calo_tree->Branch("allhit_width_v", "std::vector<float>", &_allhit_width_v);
  _calo_tree->Branch("allhit_width_y", "std::vector<float>", &_allhit_width_y);
  _calo_tree->Branch("allhit_time_u", "std::vector<float>", &_allhit_time_u);
  _calo_tree->Branch("allhit_time_v", "std::vector<float>", &_allhit_time_v);
  _calo_tree->Branch("allhit_time_y", "std::vector<float>", &_allhit_time_y);

  _calo_tree->Branch("trkhit_charge_u", "std::vector<float>", &_trkhit_charge_u);
  _calo_tree->Branch("trkhit_charge_v", "std::vector<float>", &_trkhit_charge_v);
  _calo_tree->Branch("trkhit_charge_y", "std::vector<float>", &_trkhit_charge_y);
  _calo_tree->Branch("trkhit_channel_u", "std::vector<unsigned>", &_trkhit_channel_u);
  _calo_tree->Branch("trkhit_channel_v", "std::vector<unsigned>", &_trkhit_channel_v);
  _calo_tree->Branch("trkhit_channel_y", "std::vector<unsigned>", &_trkhit_channel_y);
  _calo_tree->Branch("trkhit_multiplicity_u", "std::vector<int>", &_trkhit_multiplicity_u); 
  _calo_tree->Branch("trkhit_multiplicity_v", "std::vector<int>", &_trkhit_multiplicity_v); 
  _calo_tree->Branch("trkhit_multiplicity_y", "std::vector<int>", &_trkhit_multiplicity_y); 
  _calo_tree->Branch("trkhit_width_u", "std::vector<float>", &_trkhit_width_u);
  _calo_tree->Branch("trkhit_width_v", "std::vector<float>", &_trkhit_width_v);
  _calo_tree->Branch("trkhit_width_y", "std::vector<float>", &_trkhit_width_y);
  _calo_tree->Branch("trkhit_time_u", "std::vector<float>", &_trkhit_time_u);
  _calo_tree->Branch("trkhit_time_v", "std::vector<float>", &_trkhit_time_v);
  _calo_tree->Branch("trkhit_time_y", "std::vector<float>", &_trkhit_time_y);

  _calo_tree->Branch("calohit_charge_u", "std::vector<float>", &_calohit_charge_u);
  _calo_tree->Branch("calohit_charge_v", "std::vector<float>", &_calohit_charge_v);
  _calo_tree->Branch("calohit_charge_y", "std::vector<float>", &_calohit_charge_y);
  _calo_tree->Branch("calohit_channel_u", "std::vector<unsigned>", &_calohit_channel_u);
  _calo_tree->Branch("calohit_channel_v", "std::vector<unsigned>", &_calohit_channel_v);
  _calo_tree->Branch("calohit_channel_y", "std::vector<unsigned>", &_calohit_channel_y);
  _calo_tree->Branch("calohit_multiplicity_u", "std::vector<int>", &_calohit_multiplicity_u); 
  _calo_tree->Branch("calohit_multiplicity_v", "std::vector<int>", &_calohit_multiplicity_v); 
  _calo_tree->Branch("calohit_multiplicity_y", "std::vector<int>", &_calohit_multiplicity_y); 
  _calo_tree->Branch("calohit_width_u", "std::vector<float>", &_calohit_width_u);
  _calo_tree->Branch("calohit_width_v", "std::vector<float>", &_calohit_width_v);
  _calo_tree->Branch("calohit_width_y", "std::vector<float>", &_calohit_width_y);
  _calo_tree->Branch("calohit_time_u", "std::vector<float>", &_calohit_time_u);
  _calo_tree->Branch("calohit_time_v", "std::vector<float>", &_calohit_time_v);
  _calo_tree->Branch("calohit_time_y", "std::vector<float>", &_calohit_time_y);

  _calo_tree->Branch("sumhit_charge_u", "std::vector<float>", &_sumhit_charge_u);
  _calo_tree->Branch("sumhit_charge_v", "std::vector<float>", &_sumhit_charge_v);
  _calo_tree->Branch("sumhit_charge_y", "std::vector<float>", &_sumhit_charge_y);
  _calo_tree->Branch("sumhit_channel_u", "std::vector<unsigned>", &_sumhit_channel_u);
  _calo_tree->Branch("sumhit_channel_v", "std::vector<unsigned>", &_sumhit_channel_v);
  _calo_tree->Branch("sumhit_channel_y", "std::vector<unsigned>", &_sumhit_channel_y);
  _calo_tree->Branch("sumhit_multiplicity_u", "std::vector<int>", &_sumhit_multiplicity_u); 
  _calo_tree->Branch("sumhit_multiplicity_v", "std::vector<int>", &_sumhit_multiplicity_v); 
  _calo_tree->Branch("sumhit_multiplicity_y", "std::vector<int>", &_sumhit_multiplicity_y); 
  _calo_tree->Branch("sumhit_width_u", "std::vector<float>", &_sumhit_width_u);
  _calo_tree->Branch("sumhit_width_v", "std::vector<float>", &_sumhit_width_v);
  _calo_tree->Branch("sumhit_width_y", "std::vector<float>", &_sumhit_width_y);
  _calo_tree->Branch("sumhit_time_u", "std::vector<float>", &_sumhit_time_u);
  _calo_tree->Branch("sumhit_time_v", "std::vector<float>", &_sumhit_time_v);
  _calo_tree->Branch("sumhit_time_y", "std::vector<float>", &_sumhit_time_y);

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
  _calo_tree->Branch("daughters", &_daughters, "daughters/i");
  _calo_tree->Branch("shr_daughters", &_shr_daughters, "shr_daughters/i");
  _calo_tree->Branch("n_pfp", &_n_pfp, "n_pfp/i");

  // track information
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

  _calo_tree->Branch("sx", "std::vector<float>", & _sx); 
  _calo_tree->Branch("sy", "std::vector<float>", & _sy); 
  _calo_tree->Branch("sz", "std::vector<float>", & _sz); 

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
            const std::vector<art::Ptr<recob::Wire>> &allWires,
            const std::vector<art::Ptr<recob::Hit>> &allHits,
            const art::Ptr<recob::PFParticle> &pfp,
            const art::Ptr<recob::Track> &trk,
            const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
            const std::vector<art::Ptr<anab::Calorimetry>> &calos,
            const std::vector<art::Ptr<anab::ParticleID>> &pids,
            const std::vector<art::Ptr<recob::Hit>> &areaHits,
            const std::vector<art::Ptr<recob::Hit>> &trkHits,
            const std::vector<art::Ptr<recob::Hit>> &caloHits,
            const recob::MCSFitResult &muon_mcs,
            const sbn::RangeP &muon_range,
            const sbn::RangeP &proton_range,
            const bool fData,
            const bool fShrFit,
            const std::vector<std::pair<int, float>> &particleMatches,
            const art::Ptr<simb::MCParticle> &trueParticle,
            const std::array<std::map<unsigned, std::vector<const sim::IDE*>>, 3> particle_ide_map,
            const std::array<std::vector<const sim::IDE*>,3> &allIDEs,
            const std::vector<geo::BoxBoundedGeo> &activeVolumes)
            //const std::vector<searchingfornues::BtPart> btparts_v,
            //const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> &assocMCPart)
{

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // TODO: variables to set
  _generation = -1;
  _shr_daughters = 0;
  _trk_daughters = 0;
  _trk_score = -1;
  _longest = -1;

  _daughters = pfp->NumDaughters();

  for (unsigned sp_i = 0; sp_i < spacepoints.size(); sp_i++) {
    _sx.push_back(spacepoints[sp_i]->XYZ()[0]);
    _sy.push_back(spacepoints[sp_i]->XYZ()[1]);
    _sz.push_back(spacepoints[sp_i]->XYZ()[2]);
  }

  for (const auto &wire: allWires) {
    unsigned plane_id = geometry->ChannelToWire(wire->Channel()).at(0).Plane;

    if (plane_id == 0) {
      _integrated_charge_u.push_back(0);
      for (auto const &range: wire->SignalROI().get_ranges()) {
        for (float val: range) _integrated_charge_u.back() += val;
      }
    }
    else if (plane_id == 1) {
      _integrated_charge_v.push_back(0);
      for (auto const &range: wire->SignalROI().get_ranges()) {
        for (float val: range) _integrated_charge_v.back() += val;
      }
    }
    else if (plane_id == 2) {
      _integrated_charge_y.push_back(0);
      for (auto const &range: wire->SignalROI().get_ranges()) {
        for (float val: range) _integrated_charge_y.back() += val;
      }
    }
  }

  FillHits(allHits, 
    _allhit_charge_u,
    _allhit_charge_v,
    _allhit_charge_y,
    _allhit_channel_u,
    _allhit_channel_v,
    _allhit_channel_y,
    _allhit_multiplicity_u,
    _allhit_multiplicity_v,
    _allhit_multiplicity_y,
    _allhit_width_u,
    _allhit_width_v,
    _allhit_width_y,
    _allhit_time_u,
    _allhit_time_v,
    _allhit_time_y
  );

  FillHits(areaHits, 
    _areahit_charge_u,
    _areahit_charge_v,
    _areahit_charge_y,
    _areahit_channel_u,
    _areahit_channel_v,
    _areahit_channel_y,
    _areahit_multiplicity_u,
    _areahit_multiplicity_v,
    _areahit_multiplicity_y,
    _areahit_width_u,
    _areahit_width_v,
    _areahit_width_y,
    _areahit_time_u,
    _areahit_time_v,
    _areahit_time_y
  );

  FillHits(trkHits, 
    _trkhit_charge_u,
    _trkhit_charge_v,
    _trkhit_charge_y,
    _trkhit_channel_u,
    _trkhit_channel_v,
    _trkhit_channel_y,
    _trkhit_multiplicity_u,
    _trkhit_multiplicity_v,
    _trkhit_multiplicity_y,
    _trkhit_width_u,
    _trkhit_width_v,
    _trkhit_width_y,
    _trkhit_time_u,
    _trkhit_time_v,
    _trkhit_time_y
  );
    
  FillHits(caloHits, 
    _calohit_charge_u,
    _calohit_charge_v,
    _calohit_charge_y,
    _calohit_channel_u,
    _calohit_channel_v,
    _calohit_channel_y,
    _calohit_multiplicity_u,
    _calohit_multiplicity_v,
    _calohit_multiplicity_y,
    _calohit_width_u,
    _calohit_width_v,
    _calohit_width_y,
    _calohit_time_u,
    _calohit_time_v,
    _calohit_time_y
  );

  FillHits(caloHits, 
    _sumhit_charge_u,
    _sumhit_charge_v,
    _sumhit_charge_y,
    _sumhit_channel_u,
    _sumhit_channel_v,
    _sumhit_channel_y,
    _sumhit_multiplicity_u,
    _sumhit_multiplicity_v,
    _sumhit_multiplicity_y,
    _sumhit_width_u,
    _sumhit_width_v,
    _sumhit_width_y,
    _sumhit_time_u,
    _sumhit_time_v,
    _sumhit_time_y,
    false);

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
    _backtracked_q_u = 0.;
    _backtracked_q_v = 0.;
    _backtracked_q_y = 0.;
    _backtracked_qtotal_u = 0.;
    _backtracked_qtotal_v = 0.;
    _backtracked_qtotal_y = 0.;
    _backtracked_deposited_e_u = 0.;
    _backtracked_deposited_e_v = 0.;
    _backtracked_deposited_e_y = 0.;

    for (auto const &ide_pair: particle_ide_map[0]) {
      _backtracked_c_u.push_back(ide_pair.first);
      _backtracked_pc_q_u.push_back(0);
      for (const sim::IDE *ide: ide_pair.second) {
        trueParticleE += ide->energy;
        _backtracked_q_u += ide->numElectrons;
        _backtracked_deposited_e_u += ide->energy;
        _backtracked_pc_q_u.back() += ide->numElectrons;
      }
    }

    for (auto const &ide_pair: particle_ide_map[1]) {
      _backtracked_c_v.push_back(ide_pair.first);
      _backtracked_pc_q_v.push_back(0);
      for (const sim::IDE *ide: ide_pair.second) {
        trueParticleE += ide->energy;
        _backtracked_q_v += ide->numElectrons;
        _backtracked_deposited_e_v += ide->energy;
        _backtracked_pc_q_v.back() += ide->numElectrons;
      }
    }

    for (auto const &ide_pair: particle_ide_map[2]) {
      _backtracked_c_y.push_back(ide_pair.first);
      _backtracked_pc_q_y.push_back(0);
      for (const sim::IDE *ide: ide_pair.second) {
        trueParticleE += ide->energy;
        _backtracked_q_y += ide->numElectrons;
        _backtracked_deposited_e_y += ide->energy;
        _backtracked_pc_q_y.back() += ide->numElectrons;
      }
    }

    for (auto ide: allIDEs[0]) {
      _backtracked_qtotal_u += ide->numElectrons;
    }
    for (auto ide: allIDEs[1]) {
      _backtracked_qtotal_v += ide->numElectrons;
    }
    for (auto ide: allIDEs[2]) {
      _backtracked_qtotal_y += ide->numElectrons;
    }

    _backtracked_completeness = matchE / trueParticleE;

    // TODO: set overlay purity
    
   _backtracked_theta = trueParticle->Momentum().Theta();
   _backtracked_phi = trueParticle->Momentum().Phi();

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

    std::cout << "PID: " << _backtracked_pdg << " Process: " << trueParticle->EndProcess() << " is stopping: " << _backtracked_process_is_stopping << std::endl;

    std::vector<geoalgo::AABox> aa_volumes;
    _backtracked_end_in_tpc = false;
    for (const geo::BoxBoundedGeo &AV: activeVolumes) {
      if (AV.ContainsPosition(trueParticle->Position().Vect()) && AV.ContainsPosition(trueParticle->EndPosition().Vect())) {
        _backtracked_end_in_tpc = true;
        aa_volumes.emplace_back(AV.MinX(), AV.MinY(), AV.MinZ(), AV.MaxX(), AV.MaxY(), AV.MaxZ());
        break;
      }
    }

    _backtracked_length_start_to_end = (trueParticle->Position().Vect() - trueParticle->EndPosition().Vect()).Mag();
    _backtracked_length = 0.;
    for (unsigned i = 1; i < trueParticle->NumberTrajectoryPoints(); i++) {
      _backtracked_length += ContainedLength(trueParticle->Position(i-1).Vect(), trueParticle->Position(i).Vect(), aa_volumes);
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
    float min_dist = 1000;
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

    if (i_min == (size_t)-1) return;

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

void FillHits(const std::vector<art::Ptr<recob::Hit>> &hits,
              std::vector<float> &charge_u,
              std::vector<float> &charge_v,
              std::vector<float> &charge_y,
              std::vector<unsigned> &channel_u,
              std::vector<unsigned> &channel_v,
              std::vector<unsigned> &channel_y,
              std::vector<int> &multiplicity_u,
              std::vector<int> &multiplicity_v,
              std::vector<int> &multiplicity_y,
              std::vector<float> &width_u,
              std::vector<float> &width_v,
              std::vector<float> &width_y,
              std::vector<float> &time_u,
              std::vector<float> &time_v,
              std::vector<float> &time_y,
              bool use_integral) {

  detinfo::DetectorClocks const* dclock
          = lar::providerFrom<detinfo::DetectorClocksService>();

  struct HitInfo {
    float charge;
    int multiplicity;
    float width;
    int start;
    int end;
    float time;
    int nhit;
    HitInfo(): charge(0), multiplicity(0), width(0), start(0), end(0), time(0.), nhit(0) {}
  };

  std::map<unsigned, HitInfo> hit_map_u;
  std::map<unsigned, HitInfo> hit_map_v;
  std::map<unsigned, HitInfo> hit_map_y;

  for (const art::Ptr<recob::Hit> hit: hits) {
    if (hit->WireID().Plane == 0) {
      // if using summed ADC method to count charge, don't double-count hits in the same snippet
      if (!use_integral && hit_map_u[hit->Channel()].start == hit->StartTick() && hit_map_u[hit->Channel()].end == hit->EndTick()) {
        continue;
      }
      hit_map_u[hit->Channel()].charge += (use_integral) ? hit->Integral() : hit->SummedADC();
      hit_map_u[hit->Channel()].multiplicity = std::max(hit_map_u[hit->Channel()].multiplicity, (int) hit->Multiplicity());
      hit_map_u[hit->Channel()].width = std::max(hit_map_u[hit->Channel()].width, (float) (hit->EndTick() - hit->StartTick()));
      hit_map_u[hit->Channel()].start = hit->StartTick();
      hit_map_u[hit->Channel()].end = hit->EndTick();
      hit_map_u[hit->Channel()].time = ( dclock->TPCClock().TickPeriod() * (hit->StartTick() + hit->EndTick()) / 2. + hit_map_u[hit->Channel()].time * hit_map_u[hit->Channel()].nhit) / (hit_map_u[hit->Channel()].nhit + 1);
      hit_map_u[hit->Channel()].nhit ++;
    }
    else if (hit->WireID().Plane == 1) {
      // if using summed ADC method to count charge, don't double-count hits in the same snippet
      if (!use_integral && hit_map_v[hit->Channel()].start == hit->StartTick() && hit_map_v[hit->Channel()].end == hit->EndTick()) {
        continue;
      }
      hit_map_v[hit->Channel()].charge += (use_integral) ? hit->Integral() : hit->SummedADC();
      hit_map_v[hit->Channel()].multiplicity = std::max(hit_map_v[hit->Channel()].multiplicity, (int) hit->Multiplicity());
      hit_map_v[hit->Channel()].width = std::max(hit_map_v[hit->Channel()].width, (float) (hit->EndTick() - hit->StartTick()));
      hit_map_v[hit->Channel()].start = hit->StartTick();
      hit_map_v[hit->Channel()].end = hit->EndTick();
      hit_map_v[hit->Channel()].time = ( dclock->TPCClock().TickPeriod() * (hit->StartTick() + hit->EndTick()) / 2. + hit_map_v[hit->Channel()].time * hit_map_v[hit->Channel()].nhit) / (hit_map_v[hit->Channel()].nhit + 1);
      hit_map_v[hit->Channel()].nhit ++;
    }
    else {
      // if using summed ADC method to count charge, don't double-count hits in the same snippet
      if (!use_integral && hit_map_y[hit->Channel()].start == hit->StartTick() && hit_map_y[hit->Channel()].end == hit->EndTick()) {
        continue;
      }
      hit_map_y[hit->Channel()].charge += (use_integral) ? hit->Integral() : hit->SummedADC();
      hit_map_y[hit->Channel()].multiplicity = std::max(hit_map_y[hit->Channel()].multiplicity, (int) hit->Multiplicity());
      hit_map_y[hit->Channel()].width = std::max(hit_map_y[hit->Channel()].width, (float) (hit->EndTick() - hit->StartTick()));
      hit_map_y[hit->Channel()].start = hit->StartTick();
      hit_map_y[hit->Channel()].end = hit->EndTick();
      hit_map_y[hit->Channel()].time = ( dclock->TPCClock().TickPeriod() * (hit->StartTick() + hit->EndTick()) / 2. + hit_map_y[hit->Channel()].time * hit_map_y[hit->Channel()].nhit) / (hit_map_y[hit->Channel()].nhit + 1);
      hit_map_y[hit->Channel()].nhit ++;
    }
  }

  for (auto const &pair: hit_map_u) {
    channel_u.push_back(pair.first);
    charge_u.push_back(pair.second.charge);
    multiplicity_u.push_back(pair.second.multiplicity);
    width_u.push_back(pair.second.width);
    time_u.push_back(pair.second.time);
  }
  for (auto const &pair: hit_map_v) {
    channel_v.push_back(pair.first);
    charge_v.push_back(pair.second.charge);
    multiplicity_v.push_back(pair.second.multiplicity);
    width_v.push_back(pair.second.width);
    time_v.push_back(pair.second.time);
  }
  for (auto const &pair: hit_map_y) {
    channel_y.push_back(pair.first);
    charge_y.push_back(pair.second.charge);
    multiplicity_y.push_back(pair.second.multiplicity);
    width_y.push_back(pair.second.width);
    time_y.push_back(pair.second.time);
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

float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes) {
  static const geoalgo::GeoAlgo algo;

  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  // construct individual points
  geoalgo::Point_t p0(v0);
  geoalgo::Point_t p1(v1);

  // construct line segment
  geoalgo::LineSegment line(p0, p1);

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.Contain(p0) + box.Contain(p1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist)
    if (n_contained == 1) {
      auto intersections = algo.Intersection(line, box);
      // Because of floating point errors, it can sometimes happen
      // that there is 1 contained point but no "Intersections"
      // if one of the points is right on the edge
      if (intersections.size() == 0) {
        // determine which point is on the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        assert(p0_edge || p1_edge);
        // contained one is on edge -- can treat both as not contained
        //
        // In this case, no length
        if ((p0_edge && box.Contain(p0)) || (box.Contain(p1) && p1_edge))
          continue;
        // un-contaned one is on edge -- treat both as contained
        else if ((p0_edge && box.Contain(p1)) || (box.Contain(p0) && p1_edge)) {
	  length = (v1 - v0).Mag();
	  break;
        }
        else {
          assert(false); // bad
        }
      }
      // floating point errors can also falsely cause 2 intersection points
      //
      // in this case, one of the intersections must be very close to the 
      // "contained" point, so the total contained length will be about
      // the same as the distance between the two intersection points
      else if (intersections.size() == 2) {
        length += (intersections.at(0).ToTLorentzVector().Vect() - intersections.at(1).ToTLorentzVector().Vect()).Mag();
        continue;
      }
      // "Correct"/ideal case -- 1 intersection point
      else if (intersections.size() == 1) {
        // get TVector at intersection point
        TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
        length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() ); 
      }
      else assert(false); // bad
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      if (!(intersections.size() == 0 || intersections.size() == 2)) {
        // more floating point error fixes...
        //
        // figure out which points are near the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        // and which points are near the intersection
        TVector3 vint = intersections.at(0).ToTLorentzVector().Vect();

        bool p0_int = (v0 - vint).Mag() < tol;
        bool p1_int = (v1 - vint).Mag() < tol;
        // exactly one of them should produce the intersection
        assert((p0_int && p0_edge) != (p1_int && p1_edge));
        // void variables when assert-ions are turned off
        (void) p0_int; (void) p1_int;

        // both close to edge -- full length is contained
        if (p0_edge && p1_edge) {
          length += (v0 - v1).Mag();
        }
        // otherwise -- one of them is not on an edge, no length is contained
        else {}
      }
      // assert(intersections.size() == 0 || intersections.size() == 2);
      else if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;
}//ContainedLength

} // namespace analysis

DEFINE_ART_MODULE(analysis::CalorimetryAnalysis)
#endif
