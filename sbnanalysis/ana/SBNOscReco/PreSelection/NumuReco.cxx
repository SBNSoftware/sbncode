#include <list>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <queue>

#include "fhiclcpp/ParameterSet.h"

#include <TParameter.h>
#include <TH2D.h>
#include <TH1D.h>

#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "core/Event.hh"
#include "core/Experiment.hh"
#include "NumuReco.h"
#include "../SBNOsc/Utilities.h"
#include "../RecoUtils/RecoUtils.h"
#include "../RecoUtils/GeoUtil.h"
#include "../NumuReco/PrimaryTrack.h"
#include "../NumuReco/TruthMatch.h"
#include "ana/SBNOscReco/TriggerEmulator/PMTTrigger.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "canvas/Persistency/Provenance/ProductID.h"

#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// copied in RecoUtils here to not use art services

namespace ana {
  namespace SBNOsc {

const static TVector3 InvalidTVector3 = TVector3(-999, -999, -999);

numu::G4ProcessID GetG4ProcessID(const std::string &process_name) {
#define MATCH_PROCESS(name) if (process_name == #name) {return numu::name;}
#define MATCH_PROCESS_NAMED(strname, id) if (process_name == #strname) {return numu::id;}
  MATCH_PROCESS(primary)
  MATCH_PROCESS(CoupledTransportation)
  MATCH_PROCESS(FastScintillation)
  MATCH_PROCESS(Decay)
  MATCH_PROCESS(anti_neutronInelastic)
  MATCH_PROCESS(neutronInelastic)
  MATCH_PROCESS(anti_protonInelastic)
  MATCH_PROCESS(protonInelastic)
  MATCH_PROCESS(hadInelastic)
  MATCH_PROCESS_NAMED(kaon+Inelastic, kaonpInelastic)
  MATCH_PROCESS_NAMED(kaon-Inelastic, kaonmInelastic)
  MATCH_PROCESS_NAMED(kaon+Inelastic, kaonpInelastic)
  MATCH_PROCESS_NAMED(kaon-Inelastic, kaonmInelastic)
  MATCH_PROCESS_NAMED(sigma+Inelastic, sigmapInelastic)
  MATCH_PROCESS_NAMED(sigma-Inelastic, sigmamInelastic)
  MATCH_PROCESS_NAMED(pi+Inelastic, pipInelastic)
  MATCH_PROCESS_NAMED(pi-Inelastic, pimInelastic)
  MATCH_PROCESS(kaon0LInelastic)
  MATCH_PROCESS(kaon0SInelastic)
  MATCH_PROCESS(lambdaInelastic)
  MATCH_PROCESS(He3Inelastic)
  MATCH_PROCESS(ionInelastic)
  MATCH_PROCESS(xi0Inelastic)
  MATCH_PROCESS(alphaInelastic)
  MATCH_PROCESS(tInelastic)
  MATCH_PROCESS(dInelastic)
  MATCH_PROCESS(anti_neutronElastic)
  MATCH_PROCESS(neutronElastic)
  MATCH_PROCESS(anti_protonElastic)
  MATCH_PROCESS(protonElastic)
  MATCH_PROCESS(hadElastic)
  MATCH_PROCESS_NAMED(kaon+Elastic, kaonpElastic)
  MATCH_PROCESS_NAMED(kaon-Elastic, kaonmElastic)
  MATCH_PROCESS_NAMED(pi+Elastic, pipElastic)
  MATCH_PROCESS_NAMED(pi-Elastic, pimElastic)
  MATCH_PROCESS(conv)
  MATCH_PROCESS(phot)
  MATCH_PROCESS(annihil)
  MATCH_PROCESS(nCapture)
  MATCH_PROCESS(nKiller)
  MATCH_PROCESS(muMinusCaptureAtRest)
  MATCH_PROCESS(CoulombScat)
  MATCH_PROCESS(hBertiniCaptureAtRest)
  MATCH_PROCESS(hFritiofCaptureAtRest)
  MATCH_PROCESS(photonNuclear)
  MATCH_PROCESS(muonNuclear)
  MATCH_PROCESS(electronNuclear)
  MATCH_PROCESS(positronNuclear)
  std::cerr << "Error: Process name with no match (" << process_name << ")\n";
  assert(false);
  return numu::primary; // unreachable
#undef MATCH_PROCESS
#undef MATCH_PROCESS_NAMED
}

void DumpTrueStart(const gallery::Event &ev, int mcparticle_id) {
  // track.match.mcparticle_id);
  const std::vector<simb::MCParticle> &mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");
      std::cout << "Track: " << mcparticle_id;
  for (const simb::MCParticle &part: mcparticle_list) {
    if (mcparticle_id == part.TrackId()) {
      std::cout << " start T: " << part.Position().T() << " to: " << part.EndPosition().T() << std::endl;
    }
  }
  return;
}

numu::Wall GetWallCross(const geo::BoxBoundedGeo &volume, const TVector3 p0, const TVector3 p1) {
  TVector3 direction = (p1 - p0) * ( 1. / (p1 - p0).Mag());
  std::vector<TVector3> intersections = volume.GetIntersections(p0, direction);
  
  /*
  std::cout << "p0: " << p0.X() << " " << p0.Y() << " " << p0.Z() << std::endl;
  std::cout << "p1: " << p1.X() << " " << p1.Y() << " " << p1.Z() << std::endl;
  std::cout << "i0: " << intersections[0].X() << " " << intersections[0].Y() << " " << intersections[0].Z() << std::endl; 
  std::cout << "i1: " << intersections[1].X() << " " << intersections[1].Y() << " " << intersections[1].Z() << std::endl; 
  */

  assert(intersections.size() == 2);

  // get the intersection point closer to p0
  int intersection_i = ((intersections[0] - p0).Mag() < (intersections[1] - p0).Mag()) ? 0 : 1;
  

  double eps = 1e-3;
  if (abs(intersections[intersection_i].X() - volume.MinX()) < eps) {
    //std::cout << "Left\n";
    return numu::wLeft;
  }
  else if (abs(intersections[intersection_i].X() - volume.MaxX()) < eps) {
    //std::cout << "Right\n";
    return numu::wRight;
  }
  else if (abs(intersections[intersection_i].Y() - volume.MinY()) < eps) {
    //std::cout << "Bottom\n";
    return numu::wBottom;
  }
  else if (abs(intersections[intersection_i].Y() - volume.MaxY()) < eps) {
    //std::cout << "Top\n";
    return numu::wTop;
  }
  else if (abs(intersections[intersection_i].Z() - volume.MinZ()) < eps) {
    //std::cout << "Front\n";
    return numu::wFront;
  }
  else if (abs(intersections[intersection_i].Z() - volume.MaxZ()) < eps) {
    //std::cout << "Back\n";
    return numu::wBack;
  }
  else assert(false);
  //std::cout << "None\n";  

  return numu::wNone;
}

int NumuReco::GetPhotonMotherID(int mcparticle_id) {
  bool is_cosmic = _true_particles_to_truth.at(mcparticle_id)->Origin() != simb::Origin_t::kBeamNeutrino;
  if (is_cosmic) {
    unsigned truth_index = _true_particles_to_generator_info.at(mcparticle_id)->generatedParticleIndex();
    return (int) truth_index +1;
  }
  else {
    return mcparticle_id;
  }
}

sbnd::crt::CRTHit ICARUS2SBNDCrtHit(const icarus::crt::CRTHit& inp) {
  sbnd::crt::CRTHit ret;
  ret.feb_id = inp.feb_id;
  ret.pesmap = inp.pesmap;
  // convert ADC -> PE's
  // TODO: fix -- hardcoded for now as temporary hack
  unsigned n_strip = 2;
  double baseline = 63.6; // ADC
  double gain = 70; // ADC / PE
  ret.peshit = (inp.peshit - n_strip*baseline) / (gain * n_strip);
  ret.ts0_s = inp.ts0_s;
  ret.ts0_s_corr = inp.ts0_s_corr;
  ret.ts0_ns = inp.ts0_ns;
  ret.ts0_ns_corr = inp.ts0_ns_corr;
  ret.ts1_ns = inp.ts1_ns;
  ret.plane = inp.plane;
  ret.x_pos = inp.x_pos;
  ret.x_err = inp.x_err;
  ret.y_pos = inp.y_pos;
  ret.y_err = inp.y_err;
  ret.z_pos = inp.z_pos;
  ret.z_err = inp.z_err;
  ret.tagger = inp.tagger;
  return ret;
}

NumuReco::NumuReco() :
  SelectionBase(),
  _event_counter(0),
  _nu_count(0),
  _selected(new std::vector<numu::RecoInteraction>) {}

void NumuReco::Initialize(fhicl::ParameterSet* config) {
  if (config) {
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NumuReco");

    // configure the Cosmic ID algorithms
    _crt_track_matchalg = new sbnd::CRTTrackMatchAlg(fhicl::Table<sbnd::CRTTrackMatchAlg::Config>(pconfig.get<fhicl::ParameterSet>("CRTTrackMatchAlg"))(), 
      fProviderManager->GetGeometryProvider(), fProviderManager->GetDetectorPropertiesProvider());

    _crt_hit_matchalg = new sbnd::CRTT0MatchAlg(fhicl::Table<sbnd::CRTT0MatchAlg::Config>(pconfig.get<fhicl::ParameterSet>("CRTT0MatchAlg"))(),
      fProviderManager->GetGeometryProvider(), fProviderManager->GetDetectorPropertiesProvider());

    _apa_cross_cosmic_alg.reconfigure(*fProviderManager, fhicl::Table<ApaCrossCosmicIdAlg::Config>(pconfig.get<fhicl::ParameterSet>("ApaCrossCosmicIdAlg"))());

    _stopping_cosmic_alg.reconfigure(*fProviderManager, fhicl::Table<StoppingParticleCosmicIdAlg::Config>(pconfig.get<fhicl::ParameterSet>("StoppingParticleCosmicIdAlg"))());


    _config.tpc_volumes = SBNRecoUtils::TPCVolumes(fProviderManager->GetGeometryProvider());
    _config.active_volumes = SBNRecoUtils::ActiveVolumes(fProviderManager->GetGeometryProvider());

    // get the beam center
    _config.beamCenterX = pconfig.get<float>("beamCenterX", 130.);
    _config.beamCenterY = pconfig.get<float>("beamCenterY", 0.);

    _config.verbose = pconfig.get<bool>("verbose", false);

    _config.requireTrack = pconfig.get<bool>("requireTrack", false);

    _config.trackMatchContainmentCut = pconfig.get<double>("trackMatchContainmentCut", -1);

    _config.TSMode = pconfig.get<int>("TSMode", 1);
    _config.MakeOpHits = pconfig.get<bool>("MakeOpHits", false);
    
    _config.requireMatched = pconfig.get<bool>("requireMatched", false);
    _config.requireContained = pconfig.get<bool>("requireContained", false);
    _config.CRTHitTimeCorrection = pconfig.get<double>("CRTHitTimeCorrection", 0.);

    _config.BeamSpillWindow = pconfig.get<std::array<float, 2>>("BeamSpillWindow", {-25., 10.}); // default -- give 0.2us buffer on each side

    // setup weight config
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);
    _config.cosmicWeight = pconfig.get<double>("cosmicWeight", 1.0);

    // flash match method
    _config.FlashMatchMethod = pconfig.get<int>("FlashMatchMethod", 2);
    _config.flashMatchTimeDifference = pconfig.get<double>("flashMatchTimeDifference");

    _config.PMTTriggerThreshold = pconfig.get<int>("PMTTriggerThreshold");

    _config.CosmicIDAllTracks = pconfig.get<bool>("CosmicIDAllTracks", false);

    // whether to use flash matching in makin CRT decision
    _config.CRT2OPTimeWidth = pconfig.get<double>("CRT2OPTimeWidth", 0.);
    _config.CRTHitinOpHitRange = pconfig.get<bool>("CRTHitinOpHitRange", false);

    // get tag names
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
    _config.RecoSliceTag = config->get<std::string>("RecoSliceTag", "pandora");
    _config.RecoVertexTag = config->get<std::string>("RecoVertexTag", "pandora");
    _config.PFParticleTag = config->get<std::string>("PFParticleTag", "pandora");
    _config.FlashMatchTag = config->get<std::string>("FlashMatchTag ", "fmatch");
    _config.CaloTag = config->get<std::string>("CaloTag", "pandoraCalo");
    _config.PIDTag = config->get<std::string>("PIDTag", "pandoraPid");
    _config.CorsikaTag = config->get<std::string>("CorsikaTag", "cosmgen");
    _config.CRTTrackTag = config->get<std::string>("CRTTrackTag", "crttrack");
    _config.CRTHitTag = config->get<std::string>("CRTHitTag", "crthit");
    _config.OpFlashTag = config->get<std::string>("OpFlashTag", "ophit");
    _config.MCParticleTag = config->get<std::string>("MCParticleTag", "largeant");
    _config.TPCRecoTagSuffixes = config->get<std::vector<std::string>>("TPCRecoTagSuffixes", { "" });

    {
      fhicl::ParameterSet dCV = \
        pconfig.get<fhicl::ParameterSet>("containment_volume_inset");
      double dx = dCV.get<double>("x");
      double dy = dCV.get<double>("y");
      double zfront = dCV.get<double>("zfront");
      double zback = dCV.get<double>("zback");
      for (const geo::BoxBoundedGeo &geo: _config.active_volumes) {
        _config.containment_volumes.emplace_back(geo.MinX() + dx, geo.MaxX() - dx, geo.MinY() + dy, geo.MaxY() - dy, geo.MinZ() + zfront, geo.MaxZ() - zback);
      }
    }

    // initialize things to do with reconstruction
    double min_track_length = pconfig.get<double>("MinTrackLength", 10.);
    _track_momentum_calculator = new trkf::TrackMomentumCalculator(min_track_length);
    _mcs_fitter = new trkf::TrajectoryMCSFitter(pconfig.get<fhicl::ParameterSet>("MCSFitter", {}));

    _op_hit_maker = (_config.MakeOpHits) ? new opdet::opHitFinderSBND(pconfig.get<fhicl::ParameterSet>("OpHitMaker"), fProviderManager->GetDetectorClocksProvider()) :
        NULL;
  }

  // Setup histo's for root output
  fOutputFile->cd();

  sbnd::CRTGeoAlg crt_geo(fProviderManager->GetGeometryProvider(), fProviderManager->GetAuxDetGeometryProvider());
  std::vector<double> tagger_volumes = crt_geo.CRTLimits();

  // crt histograms
  _crt_histograms.Initialize("crt_all", tagger_volumes);

  // add branches
  fTree->Branch("reco_event", &_recoEvent);
  fTree->Branch("reco_vertices", &_selected);

  hello();
}


void NumuReco::Finalize() {
 // cleanup pointers to algorithms
 if (_crt_track_matchalg != NULL) delete _crt_track_matchalg;
 if (_crt_hit_matchalg != NULL) delete _crt_hit_matchalg;
 if (_op_hit_maker != NULL) delete _op_hit_maker;
 if (_track_momentum_calculator != NULL) delete _track_momentum_calculator;
 if (_mcs_fitter != NULL) delete _mcs_fitter;
 fOutputFile->cd();
 // finish histograms
 _crt_histograms.Write();
}

event::RecoInteraction NumuReco::CoreRecoInteraction(const std::vector<event::Interaction> &truth, const numu::RecoInteraction &vertex, double weight) {
  event::RecoInteraction ret;
  if (vertex.slice.match.mctruth_vertex_id >= 0) {
    ret.truth_index = vertex.slice.match.mctruth_vertex_id;
  }
  ret.reco_energy = vertex.nu_energy;
  ret.weight = weight;
  return ret;
}

bool NumuReco::ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &core_truth, std::vector<event::RecoInteraction>& reco) {
  if (_event_counter % 10 == 0) {
    std::cout << "NumuReco: Processing event " << _event_counter << " "
              << "(" << _nu_count << " neutrinos selected)"
              << std::endl;
  }
  // on the first event, store the type of MC we are processing
  if (_event_counter == 0) {
    // Get neutrino truth
    gallery::Handle<std::vector<simb::MCTruth>> mctruth_handle;
    bool has_mctruth = ev.getByLabel(fTruthTag, mctruth_handle);

    // also cosmics
    gallery::Handle<std::vector<simb::MCTruth>> cosmics;
    bool has_cosmics = ev.getByLabel(_config.CorsikaTag, cosmics);
    
    // cosmic + neutrino -- overlay 
    if (has_mctruth && has_cosmics) {
      fType = numu::fOverlay;
    } 
    else if (has_cosmics) {
      fType = numu::fIntimeCosmic;
    }
    else {
      fType = numu::fUnknown;
    }
    TParameter<int> mc_type("MCType", (int)fType);
    fOutputFile->cd();
    mc_type.Write();
  }


  std::cout << "New Event! " << _event_counter << std::endl;
  std::cout << "ART Event no: " << ev.eventAuxiliary().event() << std::endl;
  // clear out old containers
  _selected->clear();

  bool selected = false;

  _event_counter++;
  CollectTruthInformation(ev);

  std::map<size_t, numu::TrueParticle> particles = MCParticleInfos();

  // get the event info
  _recoEvent = Reconstruct(ev, core_truth);

  // complete the information
  _recoEvent.particles = std::move(particles);

  // set the file type
  _recoEvent.type = fType;

  // save the information
  for (unsigned i = 0; i < _recoEvent.reco.size(); i++) {
    const numu::RecoInteraction &vertex = _recoEvent.reco[i];
    // compute the weight for this interaction
    double weight = 1.;
    weight *= _config.constantWeight;
    // TODO: what about cosmics?
    // if (vertex.slice.match.mctruth_vertex_id >= 0) {
    //   for (auto const &key: _config.uniformWeights) {
    //     weight *= core_truth[vertex.slice.match.mctruth_vertex_id].weightmap.at(key)[0];
    //  }
    // }
    if (vertex.slice.match.mode == numu::mCosmic) {
      weight *= _config.cosmicWeight;
    }

    // store reco interaction
    _selected->push_back(vertex);
    // store sbncode reco interaction
    reco.push_back(CoreRecoInteraction(core_truth, vertex, weight));
    selected = true;
    _nu_count++;
  }

  return true;
}

// get information associated with true particle
numu::TrueParticle NumuReco::MCParticleInfo(const simb::MCParticle &particle) {
  numu::TrueParticle ret;

  // default values
  ret.length = 0.;
  ret.crosses_tpc = false;
  ret.wall_enter = numu::wNone;
  ret.wall_exit = numu::wNone;
  ret.deposited_energy = 0.;

  // If the interaction is outside the active volume, then g4 won't generate positions for the particle.
  // So size == 0 => outside FV
  //
  // If size != 0, then we have to check volumes
  ret.contained_in_cryo = particle.NumberTrajectoryPoints() > 0;
  ret.contained_in_tpc = particle.NumberTrajectoryPoints() > 0;
  ret.is_contained = particle.NumberTrajectoryPoints() > 0;

  // get the point at which the Trajectory enters the cryostat (which may be non-zero for cosmics)
  int entry_point = -1;

  int cryostat_index = -1;
  int tpc_index = -1;

  for (unsigned j = 0; j < particle.NumberTrajectoryPoints(); j++) {
    for (unsigned i = 0;  i < _config.active_volumes.size(); i++) {
      if (_config.active_volumes[i].ContainsPosition(particle.Position(j).Vect())) {
        entry_point = j;
        cryostat_index = i;
        break;
      }
    }
    if (entry_point != -1) break;
  }

  // find the entering wall
  if (entry_point > 0) {
    ret.wall_enter = GetWallCross(_config.active_volumes[cryostat_index], particle.Position(entry_point).Vect(), particle.Position(entry_point-1).Vect());
  }

  std::vector<geo::BoxBoundedGeo> volumes;
  if (entry_point >= 0) {
    volumes = _config.tpc_volumes[cryostat_index];
    // setup the initial TPC index
    for (int i = 0; i < volumes.size(); i++) {
      if (volumes[i].ContainsPosition(particle.Position(entry_point).Vect())) {
        tpc_index = i;
        ret.contained_in_tpc = entry_point == 0;
        break;
      }
    }
    ret.is_contained = entry_point == 0;
    ret.contained_in_cryo = entry_point == 0;
  }
  // if we couldn't find a volume for the intial point, set not contained
  else {
    ret.contained_in_cryo = false;
    ret.is_contained = false;
  }
  if (tpc_index < 0) {
    ret.contained_in_tpc = false;
  }

  // setup aa volumes too for length calc
  std::vector<geoalgo::AABox> aa_volumes;
  for (auto const &v: volumes) {
    aa_volumes.emplace_back(v.MinX(), v.MinY(), v.MinZ(), v.MaxX(), v.MaxY(), v.MaxZ());
  } 

  int exit_point = -1;

  // Get the length and determine if any point leaves the active volume
  //
  // Use every trajectory point if possible
  if (entry_point >= 0) {
    // particle trajectory
    const simb::MCTrajectory &trajectory = particle.Trajectory();
    TVector3 pos = trajectory.Position(entry_point).Vect();
    for (int i = entry_point+1; i < particle.NumberTrajectoryPoints(); i++) {
      TVector3 this_point = trajectory.Position(i).Vect();
      // get the exit point
      // update if particle is contained
      if (ret.contained_in_cryo) {
        ret.contained_in_cryo = _config.active_volumes[cryostat_index].ContainsPosition(this_point);
      }
      // check if particle has crossed TPC
      if (!ret.crosses_tpc) {
        for (int j = 0; j < volumes.size(); j++) {
          if (volumes[j].ContainsPosition(this_point) && j != tpc_index) {
            ret.crosses_tpc = true;
            break;
          }
        }
      }
      // check if particle has left tpc
      if (ret.contained_in_tpc) {
        ret.contained_in_tpc = volumes[tpc_index].ContainsPosition(this_point);
      }

      if (ret.is_contained) {
        ret.is_contained = _config.containment_volumes[cryostat_index].ContainsPosition(this_point);
      }
      
      // update length
      ret.length += containedLength(this_point, pos, aa_volumes);

      if (!_config.active_volumes[cryostat_index].ContainsPosition(this_point) && _config.active_volumes[cryostat_index].ContainsPosition(pos)) {
        exit_point = i-1;
      }

      // update energy
      if (InActive(this_point) && InActive(pos)) {
        ret.deposited_energy += trajectory.Momentum(i-1).E() - trajectory.Momentum(i).E();
      }
      pos = trajectory.Position(i).Vect();
    }
  }
  if (exit_point < 0 && entry_point >= 0) {
    exit_point = particle.NumberTrajectoryPoints() - 1; 
  }
  if (exit_point < particle.NumberTrajectoryPoints() - 1) {
    ret.wall_exit = GetWallCross(_config.active_volumes[cryostat_index], particle.Position(exit_point).Vect(), particle.Position(exit_point+1).Vect()); 
  }

  // get the generation mode
  ret.is_cosmic = _true_particles_to_truth.at(particle.TrackId())->Origin() != simb::kBeamNeutrino;

  //double costh = (entry_point >= 0 && particle.Momentum(entry_point).Vect().Mag() > 1e-4) ? particle.Pz(entry_point) / particle.Momentum(entry_point).Vect().Mag(): -999.;
  //double kinetic_energy =  (entry_point >= 0) ? particle.E(entry_point) /* already in GeV*/ - PDGMass(pdgid) / 1000. /* MeV -> GeV */: -999.;

  // other truth information
  ret.pdgid = particle.PdgCode();

  ret.start = (entry_point >= 0) ? particle.Position(entry_point).Vect(): TVector3(-9999, -9999, -9999);
  ret.start_time = (entry_point >= 0) ? particle.Position(entry_point).T() / 1000. /* ns-> us*/: -9999;
  ret.end = (exit_point >= 0) ? particle.Position(exit_point).Vect(): TVector3(-9999, -9999, -9999);
  ret.end_time = (exit_point >= 0) ? particle.Position(exit_point).T() / 1000. /* ns -> us */ : -9999;
  
  ret.start_momentum = (entry_point >= 0) ? particle.Momentum(entry_point).Vect() : TVector3(-9999, -9999, -9999);
  ret.start_energy = (entry_point >= 0) ? particle.Momentum(entry_point).E() : -9999.;
  ret.end_momentum = (exit_point >= 0) ? particle.Momentum(exit_point).Vect() : TVector3(-9999, -9999, -9999);
  ret.end_energy = (exit_point >= 0) ? particle.Momentum(exit_point).E() : -9999.;

  ret.start_process = GetG4ProcessID(particle.Process());
  ret.end_process = GetG4ProcessID(particle.EndProcess());

  ret.ID = particle.TrackId();

  return ret;
}

std::map<size_t, numu::TrueParticle> NumuReco::MCParticleInfos() {
  std::map<size_t, numu::TrueParticle> ret;
  for (unsigned i = 0; i < _true_particles.size(); i++) {
    ret[_true_particles[i]->TrackId()] = MCParticleInfo(*_true_particles[i]);
  }
  return ret;
}

double RecoTrackLength(const art::Ptr<recob::Track> &track) {
  if (track->CountValidPoints() == 0) return 0.;
  double dist = 0.;
  geo::Point_t first = track->Start();
  for (size_t i = 1; i < track->CountValidPoints(); i++) {
    geo::Point_t second = track->LocationAtPoint(i);
    dist += sqrt((second - first).Mag2());
    first = second;
  }
  return dist;
}

std::array<bool, 4> NumuReco::RecoTrackTopology(const art::Ptr<recob::Track> &track) {
  // returned info
  bool contained_in_cryo = true;
  bool contained_in_tpc = true;
  bool crosses_tpc = false; 

  bool is_contained = true;

  // start point
  geo::Point_t start = track->Start();
  // get the active volume that the start position is in
  int cryostat_index = -1;
  int tpc_index = -1;
  int containment_index = -1;
  for (int i = 0; i < _config.containment_volumes.size(); i++) {
    if (_config.containment_volumes[i].ContainsPosition(start)) {
      containment_index = i;
    }
    break;
  }
  for (int i = 0; i < _config.active_volumes.size(); i++) {
    if (_config.active_volumes[i].ContainsPosition(start)) {
      cryostat_index = i;
      break;
    }
  }
  if (containment_index < 0) {
    is_contained = false;
  }
  std::vector<geo::BoxBoundedGeo> volumes;
  if (cryostat_index >= 0) {
    volumes = _config.tpc_volumes[cryostat_index];
    for (int i = 0; i < volumes.size(); i++) {
      if (volumes[i].ContainsPosition(start)) {
        tpc_index = i;
        break;
      }
    }
  }
  else {
    contained_in_cryo = false;
  }
  if (tpc_index < 0) {
    contained_in_tpc = false;
  }

  // now check for all track points
  for (int i = 1; i < track->CountValidPoints(); i++) {
    geo::Point_t this_point = track->LocationAtPoint(i);
    if (is_contained) {
      is_contained = _config.containment_volumes[containment_index].ContainsPosition(this_point);
    }
    if (contained_in_cryo) {
      contained_in_cryo = _config.active_volumes[cryostat_index].ContainsPosition(this_point);
    }
    if (contained_in_cryo && !crosses_tpc) {
      for (int j = 0; j < volumes.size(); j++) {
        if (volumes[j].ContainsPosition(this_point) && j != tpc_index) {
          crosses_tpc = true;
          break;
        }
      }
    }
    if (contained_in_tpc) {
      contained_in_tpc = volumes[tpc_index].ContainsPosition(this_point);
    }
  }

  return {contained_in_cryo, contained_in_tpc, crosses_tpc, is_contained};
}


std::map<size_t, numu::RecoTrack> NumuReco::RecoTrackInfo() {
  std::map<size_t, numu::RecoTrack> ret;
  for (unsigned pfp_track_index = 0; pfp_track_index < _tpc_tracks.size(); pfp_track_index++) {
    const art::Ptr<recob::Track> &track = _tpc_tracks[pfp_track_index];

    // information to be saved
    numu::RecoTrack this_track;

    // Use the particle ID as a global ID
    this_track.ID = _tpc_tracks_to_particle_index.at(pfp_track_index);

    // track length
    this_track.length = track->Length();

    // get the associated PID and Calo
    assert(_tpc_tracks_to_pid.at(pfp_track_index).size() == 3); //one per plane
    
    // sum up all the pid scores weighted by n dof
    double chi2_proton = 0.;
    double chi2_kaon = 0.;
    double chi2_muon = 0.;
    double chi2_pion = 0.;
    int n_dof = 0;
    int particle_pdg = 0;
    double min_chi2 = 0.;
    for (int i =0; i < 3; i++) {
      // invalid plane means invalid calorimetry
      if (!_tpc_tracks_to_pid.at(pfp_track_index).at(i)->PlaneID()) continue;
      const art::Ptr<anab::ParticleID> &particle_id = _tpc_tracks_to_pid.at(pfp_track_index).at(i);
      // only use particle ID on collection plane
      if (fProviderManager->GetGeometryProvider()->SignalType(particle_id->PlaneID()) == geo::kCollection) {
        n_dof += particle_id->Ndf();
        chi2_proton += particle_id->Chi2Proton();
        chi2_kaon += particle_id->Chi2Kaon();
        chi2_pion += particle_id->Chi2Pion();
        chi2_muon += particle_id->Chi2Muon();
      }
    }
    if (n_dof > 0) {
      // min chi2 is PID
      std::vector<double> chi2s {chi2_proton, chi2_muon, chi2_kaon, chi2_pion};
      int min_ind = std::distance(chi2s.begin(), std::min_element(chi2s.begin(), chi2s.end()));
      min_chi2 = *std::min_element(chi2s.begin(), chi2s.end());
      if (min_ind == 0) {
        particle_pdg = 2212;
      }
      else if (min_ind == 1) {
        particle_pdg = 13;
      }
      else if (min_ind == 2) {
        particle_pdg = 312;
      }
      else if (min_ind == 3) {
        particle_pdg = 211;
      }
      else {
        assert(false);
      }
    }
    else {
      // No particle ID was provided -- set things to nonsense
      chi2_proton = -1;
      chi2_kaon = -1;
      chi2_muon = -1;
      chi2_pion = -1;
      min_chi2 = -1.5;
    }
    this_track.pdgid = particle_pdg;
    this_track.chi2_proton = chi2_proton;
    this_track.chi2_kaon = chi2_kaon;
    this_track.chi2_pion = chi2_pion;
    this_track.chi2_muon = chi2_muon;
    this_track.min_chi2 = min_chi2;
    this_track.pid_n_dof = n_dof;

    // TODO: add in calo stuff??
    assert(_tpc_tracks_to_calo.at(pfp_track_index).size() == 3);

    // calculator only has inputs for protons and muons
    this_track.range_momentum_proton = _track_momentum_calculator->GetTrackMomentum(this_track.length, 2212);
    this_track.range_momentum_muon = _track_momentum_calculator->GetTrackMomentum(this_track.length, 13);

    recob::MCSFitResult mcs_fit_muon = _mcs_fitter->fitMcs(*track, 13);
    this_track.mcs_muon.fwd_mcs_momentum = mcs_fit_muon.fwdMomentum();
    this_track.mcs_muon.fwd_mcs_momentum_err = mcs_fit_muon.fwdMomUncertainty();
    this_track.mcs_muon.bwd_mcs_momentum = mcs_fit_muon.bwdMomentum();
    this_track.mcs_muon.bwd_mcs_momentum_err = mcs_fit_muon.bwdMomUncertainty();

    recob::MCSFitResult mcs_fit_pion = _mcs_fitter->fitMcs(*track, 211);
    this_track.mcs_pion.fwd_mcs_momentum = mcs_fit_pion.fwdMomentum();
    this_track.mcs_pion.fwd_mcs_momentum_err = mcs_fit_pion.fwdMomUncertainty();
    this_track.mcs_pion.bwd_mcs_momentum = mcs_fit_pion.bwdMomentum();
    this_track.mcs_pion.bwd_mcs_momentum_err = mcs_fit_pion.bwdMomUncertainty();

    recob::MCSFitResult mcs_fit_proton = _mcs_fitter->fitMcs(*track, 2212);
    this_track.mcs_proton.fwd_mcs_momentum = mcs_fit_proton.fwdMomentum();
    this_track.mcs_proton.fwd_mcs_momentum_err = mcs_fit_proton.fwdMomUncertainty();
    this_track.mcs_proton.bwd_mcs_momentum = mcs_fit_proton.bwdMomentum();
    this_track.mcs_proton.bwd_mcs_momentum_err = mcs_fit_proton.bwdMomUncertainty();

    recob::MCSFitResult mcs_fit_kaon = _mcs_fitter->fitMcs(*track, 321);
    this_track.mcs_kaon.fwd_mcs_momentum = mcs_fit_kaon.fwdMomentum();
    this_track.mcs_kaon.fwd_mcs_momentum_err = mcs_fit_kaon.fwdMomUncertainty();
    this_track.mcs_kaon.bwd_mcs_momentum = mcs_fit_kaon.bwdMomentum();
    this_track.mcs_kaon.bwd_mcs_momentum_err = mcs_fit_kaon.bwdMomUncertainty();

    this_track.costh = track->StartDirection().Z() / sqrt( track->StartDirection().Mag2() );  

    // get track topology
    std::array<bool, 4> topology = RecoTrackTopology(track);
    this_track.crosses_tpc = topology[2];

    this_track.start = TVector3(track->Start().X(), track->Start().Y(), track->Start().Z());
    this_track.end = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    // do truth matching
    this_track.match = MatchTrack2Truth(pfp_track_index);

    // if configured, apply Cosmic ID to all tracks
    if (_config.CosmicIDAllTracks) {
      ApplyCosmicID(this_track);
    }

    ret[this_track.ID] = this_track;
  }
  return ret;
}

void NumuReco::ApplyCosmicID(numu::RecoTrack &track) {
  unsigned track_id = _tpc_particles_to_track_index.at(track.ID);
  const art::Ptr<recob::Track> &pfp_track = _tpc_tracks[track_id];

  // do CRT matching
  track.crt_match = CRTMatching(track, *pfp_track, _tpc_tracks_to_hits.at(track_id));
  
  // Cosmic ID: 
  //
  // See if start or end looks like stopping
  track.stopping_chisq_start = _stopping_cosmic_alg.StoppingChiSq(pfp_track->Vertex(), _tpc_tracks_to_calo.at(track_id));
  track.stopping_chisq_finish = _stopping_cosmic_alg.StoppingChiSq(pfp_track->End(), _tpc_tracks_to_calo.at(track_id));
}


std::vector<size_t> NumuReco::RecoSliceTracks(
    const std::map<size_t, numu::RecoTrack> &tracks,
    const std::map<size_t, numu::RecoParticle> &particles) { 

  std::vector<size_t> ret;

  for (const auto &particle_pair: particles) {
    if (tracks.count(particle_pair.first)) {
      ret.push_back(particle_pair.first);
    }
  }

  return ret;
}

std::vector<numu::RecoParticle> NumuReco::RecoParticleInfo() {
  // ret
  std::vector<numu::RecoParticle> ret;

  // iterate over all of the pfparticles
  for (size_t i = 0; i < _tpc_particles.size(); i++) {
    // new reco particle
    numu::RecoParticle this_particle;
    this_particle.ID = i;

    // get the PFParticle
    const recob::PFParticle& this_pfp = *_tpc_particles.at(i);

    // assign pandora PDG
    this_particle.pandora_pid = this_pfp.PdgCode();

    // get its daughters in the partcle "flow"
    this_particle.daughters.assign(_tpc_particles_to_daughters[i].begin(), _tpc_particles_to_daughters[i].end());
    // get the metadata
    const larpandoraobj::PFParticleMetadata& this_metadata = *_tpc_particles_to_metadata.at(i);
    // and the properties dict
    auto const &properties = this_metadata.GetPropertiesMap();
    // get the reco vertices
    for (const art::Ptr<recob::Vertex> vert: _tpc_particles_to_vertex.at(i)) {
      TVector3 vect(vert->position().X(), vert->position().Y(), vert->position().Z());
      this_particle.vertices.push_back(vect);
    }

    // access pandora special values
    if (properties.count("IsClearCosmic")) {
      this_particle.p_is_clear_cosmic = properties.at("IsClearCosmic");
    }
    else {
      this_particle.p_is_clear_cosmic = false;
    }
    if (properties.count("NuScore")) {
      this_particle.p_nu_score = properties.at("NuScore");
    }
    else {
      this_particle.p_nu_score = -1;
    }
    if (properties.count("IsNeutrino")) {
      this_particle.p_is_neutrino = properties.at("IsNeutrino");
    }
    else {
      this_particle.p_is_neutrino = false;
    }

    ret.push_back(std::move(this_particle));
  }

  return ret;
}


bool NumuReco::HasPrimaryTrack(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoSlice &slice) {
  if (slice.primary_index < 0) return false;

  std::cout << "Checking slice particle: " << slice.primary_index << std::endl;
  for (const auto &part_pair: slice.particles) {
    std::cout << "ID: " << part_pair.first << "\n";
  }
  const numu::RecoParticle &neutrino = slice.particles.at(slice.primary_index);
  for (size_t pfp_index: neutrino.daughters) {
    std::cout << "Neutrino daughter ID: " << pfp_index << std::endl;
    const numu::RecoParticle &daughter = slice.particles.at(pfp_index);
    if (tracks.count(daughter.ID)) {
      if (tracks.at(daughter.ID).length > 1e-4) {
        return true;
      }
    }
  }

  return false;
}

bool NumuReco::SelectSlice(const numu::RecoSlice &slice) {
  return slice.primary_index >= 0 && 
         slice.particles.at(slice.primary_index).p_is_neutrino && 
         abs(slice.particles.at(slice.primary_index).pandora_pid) == 14;
}

numu::TrackTruthMatch NumuReco::MatchTrack2Truth(size_t pfp_track_id) {
  // get the service
  const cheat::ParticleInventory *inventory_service = fProviderManager->GetParticleInventoryProvider();

  // get its hits
  std::vector<art::Ptr<recob::Hit>> hits = _tpc_tracks_to_hits.at(pfp_track_id);

  // this id is the same as the mcparticle ID as long as we got it from geant4
  // int id = SBNRecoUtils::TrueParticleIDFromTotalRecoHits(*fProviderManager, hits, false);
  int mcp_track_id = SBNRecoUtils::TrueParticleIDFromTotalTrueEnergy(*fProviderManager, hits, true);


  int mcparticle_index = -1;
  for (int i = 0; i < _true_particles.size(); i++) {
    if (_true_particles[i]->TrackId() == mcp_track_id) {
      mcparticle_index = i;
      break;
    }
  }

  // number returned to mean "NULL"
  // no match
  if (mcparticle_index == -1) return numu::TrackTruthMatch();

  // We got a match! Try to identify it to an origin
  art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(mcp_track_id);
  // and calculate the completion
  double completion = SBNRecoUtils::TrackCompletion(*fProviderManager, mcp_track_id, hits);

  numu::TrackTruthMatch ret;
  // ret.mctruth = truth.get();
  ret.has_match = true;
  ret.mctruth_has_neutrino = truth->NeutrinoSet();
  ret.mctruth_vertex = ret.mctruth_has_neutrino ? truth->GetNeutrino().Nu().Position().Vect() : TVector3(-999, -999, -999);
  ret.mctruth_origin = truth->Origin();

  // TODO: fix -- in time cosmics will not have the origin set correctly. Fix and set this to 2.
  if (ret.mctruth_origin == simb::kUnknown) {
    ret.mctruth_origin = simb::kCosmicRay; 
  }

  ret.mctruth_ccnc = ret.mctruth_has_neutrino ? truth->GetNeutrino().CCNC() : -1;
  ret.mcparticle_id = mcp_track_id;
  ret.completion = completion;
  ret.purity = SBNRecoUtils::TrackPurity(*fProviderManager, mcp_track_id, hits);
  ret.match_pdg = _true_particles[mcparticle_index]->PdgCode(); 
  ret.is_primary = _true_particles[mcparticle_index]->Process() == "primary"; 
  return ret;
}

std::vector<numu::RecoSlice> NumuReco::RecoSliceInfo(
    std::map<size_t, numu::RecoTrack> &reco_tracks,
    const std::vector<numu::RecoParticle> &particles) {

  std::vector<numu::RecoSlice> ret;
  // store all the particles in the slice
  for (size_t i = 0; i < _tpc_slices.size(); i++) {
    const recob::Slice &slice = *_tpc_slices[i];
    numu::RecoSlice slice_ret;
    slice_ret.primary_index = -1;
    bool primary_particle_set = false;
    // get the particles
    const std::vector<art::Ptr<recob::PFParticle>> &pfp_particles = _tpc_slices_to_particles.at(i);
    // find the primary one and store all of them
    // Find the priamry particle
    for (unsigned j = 0; j < pfp_particles.size(); j++) {
      const art::Ptr<recob::PFParticle> pfp_part = pfp_particles[j];
      if (pfp_part->IsPrimary()) {
        assert(!primary_particle_set);
        primary_particle_set = true;
        slice_ret.primary_index = _tpc_slices_to_particle_index[i][j];
      }
    }
    // match the RecoParticle's to the once that match this slice
    for (const numu::RecoParticle &part: particles) {
      if (std::find(_tpc_slices_to_particle_index[i].begin(), _tpc_slices_to_particle_index[i].end(), (unsigned)part.ID) != _tpc_slices_to_particle_index[i].end()) {
        slice_ret.particles[part.ID] = part;
      }
    } 
    // throw away bad slices
    if (!SelectSlice(slice_ret)) continue;
    
    // now get information from particles which are tracks
    slice_ret.tracks = RecoSliceTracks(reco_tracks, slice_ret.particles);

    std::cout << "Primary index: " << slice_ret.primary_index << std::endl;
    std::cout << "Is primary: " << _tpc_particles[slice_ret.primary_index]->IsPrimary() << std::endl;
    if (_tpc_particles_to_flashT0.size() && _tpc_particles_to_flashT0.at(slice_ret.primary_index).size()) {
      const anab::T0 &fmatch = *_tpc_particles_to_flashT0.at(slice_ret.primary_index).at(0);
      slice_ret.flash_match.present = true;
      slice_ret.flash_match.time = fmatch.Time(); 
      slice_ret.flash_match.score = fmatch.TriggerConfidence();
      slice_ret.flash_match.pe = fmatch.TriggerType();
      std::cout << "Match time: " << slice_ret.flash_match.time << std::endl;
      std::cout << "Match score: " << slice_ret.flash_match.score << std::endl;
      std::cout << "Match PE: " << slice_ret.flash_match.pe << std::endl;
    }
    else {
      slice_ret.flash_match.present = false;
      std::cout << "No match :(\n";
    }

    // throw away slices which do not have a primary track candidate
    if (!HasPrimaryTrack(reco_tracks, slice_ret)) continue;
    
    // if we didn't do the cosmic ID for all tracks, do it for all priamry track candidates
    if (!_config.CosmicIDAllTracks) {
      const numu::RecoParticle &neutrino = slice_ret.particles.at(slice_ret.primary_index);
      for (size_t ind: neutrino.daughters) {
        const numu::RecoParticle &daughter = slice_ret.particles.at(ind);
        if (reco_tracks.count(daughter.ID)) {
          ApplyCosmicID(reco_tracks.at(daughter.ID));
        }
      }
    }

    ret.push_back(std::move(slice_ret));
  }

  return ret;
}

void NumuReco::CollectTruthInformation(const gallery::Event &ev) {
  _true_particles.clear();
  _true_particles_to_truth.clear();
  _true_particles_to_generator_info.clear();

  // mc particles
  auto const &mcparticles = ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
  art::fill_ptr_vector(_true_particles, mcparticles);

  // MC particle (G4) to MCTruth (generator -- corsika or genie)
  art::FindManyP<simb::MCTruth, sim::GeneratedParticleInfo> particles_to_truth(mcparticles, ev, "largeant");
  for (unsigned i = 0; i < mcparticles->size(); i++) {
    _true_particles_to_truth[mcparticles->at(i).TrackId()] = particles_to_truth.at(i).at(0);
    _true_particles_to_generator_info[mcparticles->at(i).TrackId()] = particles_to_truth.data(i).at(0);
  }
}

void NumuReco::CollectTPCInformation(const gallery::Event &ev) {
  // clear out old stuff
  _tpc_slices.clear();
  _tpc_tracks.clear();
  _tpc_particles.clear();
  _tpc_slices_to_particles.clear();
  _tpc_slices_to_particle_index.clear();
  _tpc_tracks_to_particles.clear();
  _tpc_tracks_to_particle_index.clear();
  _tpc_particles_to_track_index.clear();
  _tpc_tracks_to_calo.clear();
  _tpc_tracks_to_pid.clear();
  _tpc_tracks_to_hits.clear();
  _tpc_particles_to_T0.clear();
  _tpc_particles_to_vertex.clear();
  _tpc_particles_to_daughters.clear();
  _tpc_particles_to_metadata.clear();
  _tpc_particles_to_flashT0.clear();

  for (const std::string &suffix: _config.TPCRecoTagSuffixes) {
    // Offset to convert a pandora particle ID to a NumuReco particle ID
    //
    // The different ID's are required because in ICARUS TPC reconstruction
    // is done per-cryostat. Thus, the pandora particle ID's may not be unique
    unsigned particle_id_offset = _tpc_particles.size();

    // get the slices
    const auto &slice_handle = ev.getValidHandle<std::vector<recob::Slice>>(_config.RecoSliceTag + suffix);
    art::fill_ptr_vector(_tpc_slices, slice_handle);

    // get the particles
    const auto &particle_handle = ev.getValidHandle<std::vector<recob::PFParticle>>(_config.PFParticleTag + suffix);
    art::fill_ptr_vector(_tpc_particles, particle_handle);
   
    // get the tracks
    const auto &track_handle = ev.getValidHandle<std::vector<recob::Track>>(_config.RecoTrackTag + suffix);
    art::fill_ptr_vector(_tpc_tracks, track_handle);

    // slice to particle
    art::FindManyP<recob::PFParticle> slice_to_particle(slice_handle, ev, _config.PFParticleTag + suffix); 
    unsigned i0 = _tpc_slices_to_particle_index.size();
    for (unsigned i = 0; i < slice_handle->size(); i++) {
      unsigned iset = i + i0;
      _tpc_slices_to_particles.push_back(slice_to_particle.at(i));
      _tpc_slices_to_particle_index.emplace_back();
      for (unsigned j = 0; j < slice_to_particle.at(i).size(); j++) {
        _tpc_slices_to_particle_index[iset].push_back(slice_to_particle.at(i)[j]->Self() + particle_id_offset);
        assert(_tpc_particles[_tpc_slices_to_particle_index[iset][j]] == _tpc_slices_to_particles[iset][j]);
      } 
    }
  
    // track to particle and particle to track
    art::FindManyP<recob::PFParticle> track_to_particle(track_handle, ev, _config.RecoTrackTag + suffix);
    i0 = _tpc_tracks_to_particles.size();
    for (unsigned i = 0; i < track_handle->size(); i++) {
      unsigned iset = i + i0;
      _tpc_tracks_to_particles.push_back(track_to_particle.at(i).at(0));
      _tpc_tracks_to_particle_index.push_back(particle_id_offset + _tpc_tracks_to_particles[iset]->Self());
      _tpc_particles_to_track_index[particle_id_offset + _tpc_tracks_to_particles[iset]->Self()] = iset;
      // if this isn't true something went wrong
      assert(_tpc_particles[_tpc_tracks_to_particle_index[iset]] == _tpc_tracks_to_particles[iset]);
    }

    // track to Calo
    art::FindManyP<anab::Calorimetry> track_to_calo(track_handle, ev, _config.CaloTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_calo.push_back(track_to_calo.at(i));
    }

    // track to PID
    art::FindManyP<anab::ParticleID> track_to_pid(track_handle, ev, _config.PIDTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_pid.push_back(track_to_pid.at(i));
    }

    // track to hits
    art::FindManyP<recob::Hit> track_to_hits(track_handle, ev, _config.RecoTrackTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_hits.push_back(track_to_hits.at(i));
    }

    // particle to T0 
    art::FindManyP<anab::T0> particle_to_T0(particle_handle, ev, _config.PFParticleTag + suffix);
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      _tpc_particles_to_T0.push_back(particle_to_T0.at(i));
    }

    gallery::Handle<art::Assns<recob::PFParticle,anab::T0,void>> assns_check;
    ev.getByLabel(_config.FlashMatchTag + suffix, assns_check);
    if (assns_check.isValid()) {
      // particle to flash match
      art::FindManyP<anab::T0> particle_to_flash(particle_handle, ev, _config.FlashMatchTag + suffix);
      for (unsigned i = 0; i < particle_handle->size(); i++) {
        assert(particle_to_flash.at(i).size() == 0 || particle_to_flash.at(i).size() == 1);
        _tpc_particles_to_flashT0.push_back(particle_to_flash.at(i));
      }
    }

    // particle to vertex
    art::FindManyP<recob::Vertex> particle_to_vertex(particle_handle, ev, _config.PFParticleTag + suffix);
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      assert(particle_to_vertex.at(i).size() == 0 || particle_to_vertex.at(i).size() == 1);
      _tpc_particles_to_vertex.push_back(particle_to_vertex.at(i));
    }

    // particle to daughters
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      std::vector<unsigned> daughters;
      for (unsigned d: _tpc_particles[i+particle_id_offset]->Daughters()) {
        daughters.push_back(particle_id_offset + d); 
      }
      _tpc_particles_to_daughters.push_back(std::move(daughters));
    }

    // particle to metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_metadatas(particle_handle, ev, _config.PFParticleTag + suffix);
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      _tpc_particles_to_metadata.push_back(pfp_metadatas.at(i).at(0));
    }
  }
}


void NumuReco::CollectCRTInformation(const gallery::Event &ev) {
  _crt_tracks_local.clear();
  _crt_hits_local.clear();

  // collect the CRT Tracks
  gallery::Handle<std::vector<sbnd::crt::CRTTrack>> crt_tracks_sbnd;
  bool has_crt_tracks_sbnd = ev.getByLabel(_config.CRTTrackTag, crt_tracks_sbnd);
  gallery::Handle<std::vector<icarus::crt::CRTTrack>> crt_tracks_icarus;
  // bool has_crt_tracks_icarus = ev.getByLabel(_config.CRTTrackTag, crt_tracks_icarus);
  bool has_crt_tracks_icarus = false;
  _has_crt_tracks = has_crt_tracks_sbnd || has_crt_tracks_icarus;

  if (has_crt_tracks_icarus) {
    for (unsigned i = 0; i < crt_tracks_icarus->size(); i++) {
      _crt_tracks_local.emplace_back();
      memcpy(&_crt_tracks_local[i], &crt_tracks_icarus->at(i), sizeof(sbnd::crt::CRTTrack));
    }  
  }
  _crt_tracks = (_has_crt_tracks) ?
    ( (has_crt_tracks_sbnd) ? crt_tracks_sbnd.product() : &_crt_tracks_local) : 
    NULL;

  // collect the CRT Hits
  gallery::Handle<std::vector<sbnd::crt::CRTHit>> crt_hits_sbnd;
  bool has_crt_hits_sbnd = ev.getByLabel(_config.CRTHitTag, crt_hits_sbnd);
  gallery::Handle<std::vector<icarus::crt::CRTHit>> crt_hits_icarus;
  bool has_crt_hits_icarus = ev.getByLabel(_config.CRTHitTag, crt_hits_icarus);
  _has_crt_hits = has_crt_hits_sbnd || has_crt_hits_icarus;

  // Convert ICARUS Hits to SBND Hits
  if (has_crt_hits_icarus) {
    for (unsigned i = 0; i < crt_hits_icarus->size(); i++) {
      _crt_hits_local.push_back(ICARUS2SBNDCrtHit(crt_hits_icarus->at(i)));
    }  
  }

  _crt_hits = (_has_crt_hits) ?
    ( (has_crt_hits_sbnd) ? crt_hits_sbnd.product() : &_crt_hits_local) :
    NULL;
}

numu::CRTMatch NumuReco::CRTMatching(
   const numu::RecoTrack &track,
   const recob::Track &pandora_track, 
   const std::vector<art::Ptr<recob::Hit>> &hits) {
  
  numu::CRTMatch match;

  std::pair<sbnd::crt::CRTTrack, double> closest_crt_track = { {}, -99999 };
  // try to find a match to CRT Track -- if we have one
  if (_has_crt_tracks) {
    closest_crt_track = _crt_track_matchalg->ClosestCRTTrackByAngle(pandora_track, hits, *_crt_tracks);
  }

  // if there is a track match, fill out the info
  if (closest_crt_track.second > -10000 /* Null value */) {
    match.track.present = true;
    if (_config.TSMode == 0) match.track.time = closest_crt_track.first.ts0_ns / 1000. /* ns -> us */;
    else match.track.time = closest_crt_track.first.ts1_ns / 1000. /* ns -> us */;
    match.track.angle = closest_crt_track.second;
  }
  else {
    match.track.present = false;
    match.track.time = -99999.;
    match.track.angle = -99999.;
  }

  // try to match a hit
  std::pair<sbnd::crt::CRTHit, double> hit_pair = std::pair<sbnd::crt::CRTHit, double>({sbnd::crt::CRTHit(), -1});
  if (_has_crt_hits) {
    numu::FlashMatch flash_match; // TODO: fix
    if (_config.CRTHitinOpHitRange && flash_match.present) {
      const numu::FlashMatch &match = flash_match;
      double time_width = (_config.CRT2OPTimeWidth) / 2;
      std::pair<double, double> time_range; 
      time_range.first = match.time - time_width;
      time_range.second = match.time + time_width;
      int drift_dir = sbnd::TPCGeoUtil::DriftDirectionFromHits(fProviderManager->GetGeometryProvider(), hits);
      hit_pair = _crt_hit_matchalg->ClosestCRTHit(pandora_track, time_range, *_crt_hits, drift_dir); 
    }
    else {
      std::cout << "Tryin CRT Hit match\n";
      hit_pair = _crt_hit_matchalg->ClosestCRTHit(pandora_track, hits, *_crt_hits);
    }
  }

  double distance = hit_pair.second;
  if (distance >= 0 /* matching succeeded*/) {
    match.hit_match.present = true;
    match.hit_match.distance = distance;
    match.hit = SBND2numuCRTHit(hit_pair.first);
    match.hit_match.time = match.hit.time;
    std::cout << "Match: " << match.hit.location.X() << std::endl;
  }
  else {
    match.hit_match.present = false;
    match.hit_match.distance = -99999.;
  }
 
  return match;
}

void NumuReco::CollectPMTInformation(const gallery::Event &ev) {
  // std::vector<art::Ptr<recob::OpHit>> op_hit_ptrs;
  // std::vector<recob::OpHit> op_hits;

  _op_hit_ptrs.clear();
  _op_hits_local.clear();

  return;

   // get the OpDet hits
  art::ProductID local_opdets("local_opdets");
  if (_config.MakeOpHits) {
    const std::vector<raw::OpDetWaveform> &op_waveforms = *ev.getValidHandle<std::vector<raw::OpDetWaveform>>("opdaq");
    _op_hits_local = _op_hit_maker->MakeHits(op_waveforms);
    for (unsigned i = 0; i < _op_hits_local.size(); i++) {
      _op_hit_ptrs.emplace_back(local_opdets, &_op_hits_local[i], i);
    }
  }
  else {
    gallery::Handle<std::vector<recob::OpHit>> op_hits;
    ev.getByLabel(_config.OpFlashTag, op_hits);
    art::fill_ptr_vector(_op_hit_ptrs, op_hits);
  }

}


numu::FlashMatch NumuReco::FlashMatching(
  const recob::Track &pandora_track,
  const numu::RecoTrack &track) {

  numu::FlashMatch ret;
  ret.present = false;

  // don't run if back tracker is not setup
  if (fProviderManager->GetPhotonBackTrackerProvider() == NULL) {
    return ret;
  }

  // see if we can do the match
  if (track.match.has_match) {
    // DumpTrueStart(ev, track.match.mcparticle_id);
    // get the correct index
    int photon_mother_id = GetPhotonMotherID(track.match.mcparticle_id);
    const std::vector<art::Ptr<recob::OpHit>> hit_matches = fProviderManager->GetPhotonBackTrackerProvider()->TrackIdToOpHits_Ps(photon_mother_id, _op_hit_ptrs);
    if (hit_matches.size() > 0) {
      // average the times weighted by PE
      double average = 0.;
      double match_pe = 0.;
      double min_time = 99999999.;
      for (const art::Ptr<recob::OpHit> op_hit: hit_matches) {
        average += op_hit->PE() * op_hit->PeakTime();
        match_pe += op_hit->PE();
        if (op_hit->PeakTime() < min_time) {
          min_time = op_hit->PeakTime();
        }
      }
      if (match_pe > 0) {
        double match_time = average / match_pe;
        numu::FlashMatch match;
        match.time = min_time;
        return match;
      }
    }
  }
  return ret;
}

bool NumuReco::InBeamSpill(float time) {
  return time > _config.BeamSpillWindow[0] && time < _config.BeamSpillWindow[1];
}

numu::CRTHit NumuReco::SBND2numuCRTHit(const sbnd::crt::CRTHit &hit) {
  numu::CRTHit ret;
  if (_config.TSMode == 0) ret.time = (int)hit.ts0_ns / 1000. /* ns -> us */;
  else ret.time = (int)hit.ts1_ns / 1000. /* ns -> us */;
  ret.time += _config.CRTHitTimeCorrection;
  ret.pes = hit.peshit;
  ret.location = TVector3(hit.x_pos, hit.y_pos, hit.z_pos);
  ret.uncertainty = TVector3(hit.x_err, hit.y_err, hit.z_err);
  return ret;
}

void NumuReco::FillCRTHits() {
  for (const sbnd::crt::CRTHit &hit: *_crt_hits) {
    _crt_histograms.Fill(hit);
  }
}

std::vector<numu::CRTHit> NumuReco::InTimeCRTHits() {
  std::vector<numu::CRTHit> ret;

  for (const sbnd::crt::CRTHit &hit: *_crt_hits) {
    numu::CRTHit this_hit = NumuReco::SBND2numuCRTHit(hit);
    if (InBeamSpill(this_hit.time)) {
      ret.push_back(this_hit);
    } 
  }

  return std::move(ret);
}

numu::RecoEvent NumuReco::Reconstruct(const gallery::Event &ev, const std::vector<event::Interaction> &truth) {
  // setup products
  CollectCRTInformation(ev);
  CollectPMTInformation(ev);
  CollectTPCInformation(ev);

  // colect PFParticle information
  std::vector<numu::RecoParticle> reco_particles = RecoParticleInfo();

  // collect track information
  std::map<size_t, numu::RecoTrack> reco_tracks = RecoTrackInfo();

  // collect Pandora slice information
  std::vector<numu::RecoSlice> reco_slices = RecoSliceInfo(reco_tracks, reco_particles);

  std::vector<numu::RecoInteraction> reco;
  for (unsigned reco_i = 0; reco_i < reco_slices.size(); reco_i++) {
    const numu::RecoParticle &neutrino = reco_slices[reco_i].particles.at(reco_slices[reco_i].primary_index);

    numu::RecoInteraction this_interaction;

    this_interaction.slice = reco_slices[reco_i];
    this_interaction.position = neutrino.vertices[0];
    // First guess of primary track
    this_interaction.primary_track_index = numu::SelectLongestTrack(reco_tracks, this_interaction.slice);
    assert(this_interaction.primary_track_index >= 0);

    // TODO: get the enrgy
    this_interaction.nu_energy = -1;

    // Track multiplicity
    this_interaction.multiplicity = this_interaction.slice.particles.at(this_interaction.slice.primary_index).daughters.size();

    // leave per-particle multiplicity unset until particle ID is done
    this_interaction.npion = -1;
    this_interaction.nproton = -1;
    this_interaction.nkaon = -1;

    // do initial truth matching
    this_interaction.slice.match = numu::InteractionTruthMatch(truth, reco_tracks, this_interaction);

    // store
    reco.push_back(std::move(this_interaction));
  }

  numu::RecoEvent event;
  event.reco = std::move(reco);
  event.tracks = std::move(reco_tracks);

  // collect spare detector information
  event.in_time_crt_hits = InTimeCRTHits();
  gallery::Handle<std::vector<raw::OpDetWaveform>> waveforms;
  if (ev.getByLabel("opdaq", waveforms)) {
    double tick_period = fProviderManager->GetDetectorClocksProvider()->OpticalClock().TickPeriod();
    int threshold = _config.PMTTriggerThreshold; 
    bool is_sbnd = fExperimentID == kExpSBND;
    std::pair<double, double> window;
    if (is_sbnd) window = {0., 2.}; // go out to 0.4us past end of beam gate
    else window = {1500., 1502.};
    event.flash_trigger_primitives = numu::TriggerPrimitives(*waveforms, tick_period, window, threshold, is_sbnd);
  }

  // fill spare histograms
  FillCRTHits();
  
 
  return std::move(event);
}

bool NumuReco::InActive(const TVector3 &v) const {
  for (auto const& active: _config.active_volumes) {
    if (active.ContainsPosition(v)) return true;
  }
  return false;
}


  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuReco)

