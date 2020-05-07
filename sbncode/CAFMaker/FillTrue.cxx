#include "FillTrue.h"

#include "larcorealg/GeoAlgo/GeoAlgo.h"
#include "RecoUtils/RecoUtils.h"

#include <functional>
#include <algorithm> 

// helper function declarations
caf::Wall_t GetWallCross( const geo::BoxBoundedGeo &volume, 
			  const TVector3 p0, 
			  const TVector3 p1);

caf::g4_process_ GetG4ProcessID(const std::string &name);

caf::SRTrackTruth MatchTrack2Truth(const std::vector<art::Ptr<recob::Hit>> &hits);

caf::SRTruthMatch MatchSlice2Truth(const std::vector<art::Ptr<recob::Hit>> &hits,
					  const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
					  const cheat::ParticleInventoryService &inventory_service);

float ContainedLength(const TVector3 &v0, const TVector3 &v1,
                      const std::vector<geoalgo::AABox> &boxes);

bool FRFillNumuCC(const simb::MCTruth &mctruth, 
                  const std::vector<art::Ptr<sim::MCTrack>> &mctracks, 
                  const std::vector<geo::BoxBoundedGeo> &volumes, 
                  TRandom &rand,
                  caf::SRFakeReco &fakereco);

namespace caf {

  void FillTrackTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
		      caf::SRTrack& srtrack,
		      bool allowEmpty)
  {
    // Truth matching
    srtrack.truth = MatchTrack2Truth(hits);

  }//FillTrackTruth

  //------------------------------------------------

  void FillSliceTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                      const std::vector<caf::SRTrueInteraction> &srneutrinos,
                      const cheat::ParticleInventoryService &inventory_service,
                      caf::SRSlice &srslice, caf::SRTruthBranch &srmc,
                      bool allowEmpty)
  {

    caf::SRTruthMatch tmatch;
    //    srslice.tmatch = MatchSlice2Truth(hits, neutrinos, inventory_service);
    tmatch = MatchSlice2Truth(hits, neutrinos, inventory_service);

    if (tmatch.index >= 0) {

      caf::SRTrueInteraction numatch = srneutrinos[tmatch.index];
      numatch.visEinslc         = tmatch.visEinslc;
      numatch.visEcosmic        = tmatch.visEcosmic;
      numatch.eff               = tmatch.eff;
      numatch.pur               = tmatch.pur;
      numatch.is_numucc_primary = tmatch.is_numucc_primary;

      //      srslice.truth = numatch;
      srmc.nu.push_back(numatch);
      ++srmc.nnu;
    }

    std::cout << "Slice matched to index: " << tmatch.index 
	      << " with match frac: " << tmatch.pur << std::endl;

  }//FillSliceTruth

  //------------------------------------------------
  void FillTrueNeutrino(const art::Ptr<simb::MCTruth> mctruth, 
			const art::Ptr<simb::MCFlux>  mcflux,
			std::vector<caf::SRTrueParticle> srprimaries,
			caf::SRTrueInteraction &srneutrino, size_t i) {

    srneutrino.index = i;

    if (mctruth->NeutrinoSet()) {
      // Neutrino
      const simb::MCNeutrino& nu = mctruth->GetNeutrino();
      srneutrino.isnc =   nu.CCNC()  && (nu.Mode() != simb::kWeakMix);
      srneutrino.iscc = (!nu.CCNC()) && (nu.Mode() != simb::kWeakMix);
      srneutrino.pdg = nu.Nu().PdgCode();
      srneutrino.initpdg = mcflux->fntype;
      srneutrino.targetPDG = nu.Target();
      srneutrino.genie_intcode = nu.Mode();
      srneutrino.bjorkenX = nu.X();
      srneutrino.inelasticityY = nu.Y();
      srneutrino.Q2 = nu.QSqr();
      srneutrino.w = nu.W();
      srneutrino.E = nu.Nu().EndMomentum().Energy();
      srneutrino.momentum = nu.Nu().EndMomentum().Vect();
      srneutrino.position = nu.Nu().Position().Vect();

      const simb::MCParticle& lepton = nu.Lepton();
      TLorentzVector q_labframe;
      q_labframe = nu.Nu().EndMomentum() - lepton.Momentum(0);
      srneutrino.q0_lab = q_labframe.E();
      srneutrino.modq_lab = q_labframe.P();

      srneutrino.prim = srprimaries;

    }

  }

  //------------------------------------------------

  void FillTrueG4Particle(const simb::MCParticle &particle, 
			  const std::vector<geo::BoxBoundedGeo> &active_volumes,
			  const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
                          const std::map<int, std::vector<const sim::IDE *>> &id_to_ide_map,
			  const cheat::BackTrackerService &backtracker,
			  const cheat::ParticleInventoryService &inventory_service,
			  const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                          caf::SRTrueParticle &srparticle) {

    std::vector<const sim::IDE*> empty;
    const std::vector<const sim::IDE*> &particle_ides = id_to_ide_map.count(particle.TrackId()) ? id_to_ide_map.at(particle.TrackId()) : empty;

    srparticle.length = 0.;
    srparticle.crosses_tpc = false;
    srparticle.wallin = caf::kWallNone;
    srparticle.wallout = caf::kWallNone;
    srparticle.planeVisE = 0.;
    for (auto const ide: particle_ides) {
      srparticle.planeVisE += ide->energy / 1000. /* MeV -> GeV*/;
    }

    // if no trajectory points, then assume outside AV
    srparticle.cont_tpc = particle.NumberTrajectoryPoints() > 0;
    srparticle.contained = particle.NumberTrajectoryPoints() > 0;

    // Get the entry and exit points
    int entry_point = -1;

    int cryostat_index = -1;
    int tpc_index = -1;

    for (unsigned j = 0; j < particle.NumberTrajectoryPoints(); j++) {
      for (unsigned i = 0; i < active_volumes.size(); i++) {
        if (active_volumes.at(i).ContainsPosition(particle.Position(j).Vect())) {
          entry_point = j;
          cryostat_index = i;
          break;
        }
      }
      if (entry_point != -1) break;
    }
    // get the wall
    if (entry_point > 0) {
      srparticle.wallin = GetWallCross(active_volumes.at(cryostat_index), particle.Position(entry_point).Vect(), particle.Position(entry_point-1).Vect());
    }

    int exit_point = -1;

    // now setup the cryostat the particle is in
    std::vector<geo::BoxBoundedGeo> volumes;
    if (entry_point >= 0) {
      volumes = tpc_volumes.at(cryostat_index);
      for (unsigned i = 0; i < volumes.size(); i++) {
        if (volumes[i].ContainsPosition(particle.Position(entry_point).Vect())) {
          tpc_index = i;
          srparticle.cont_tpc = entry_point == 0;
          break;
        }
      }
      srparticle.contained = entry_point == 0;
    }
    // if we couldn't find the initial point, set not contained
    else {
      srparticle.contained = false;
    }
    if (tpc_index < 0) {
      srparticle.cont_tpc = false;
    }

    // setup aa volumes too for length calc
    // Define the volume used for length calculation to be the cryostat volume in question
    std::vector<geoalgo::AABox> aa_volumes;
    if (entry_point >= 0) {
      const geo::BoxBoundedGeo &v = active_volumes.at(cryostat_index);
      aa_volumes.emplace_back(v.MinX(), v.MinY(), v.MinZ(), v.MaxX(), v.MaxY(), v.MaxZ());
    }

    // Get the length and determine if any point leaves the active volume
    //
    // Use every trajectory point if possible
    if (entry_point >= 0) {
      // particle trajectory
      const simb::MCTrajectory &trajectory = particle.Trajectory();
      TVector3 pos = trajectory.Position(entry_point).Vect();
      for (unsigned i = entry_point+1; i < particle.NumberTrajectoryPoints(); i++) {
        TVector3 this_point = trajectory.Position(i).Vect();
        // get the exit point
        // update if particle is contained
        // check if particle has crossed TPC
        if (!srparticle.crosses_tpc) {
          for (unsigned j = 0; j < volumes.size(); j++) {
            if (volumes[j].ContainsPosition(this_point) && tpc_index >= 0 && j != ((unsigned)tpc_index)) {
              srparticle.crosses_tpc = true;
              break;
            }
          }
        }
        // check if particle has left tpc
        if (srparticle.cont_tpc) {
          srparticle.cont_tpc = volumes[tpc_index].ContainsPosition(this_point);
        }

        if (srparticle.contained) {
          srparticle.contained = active_volumes.at(cryostat_index).ContainsPosition(this_point);
        }
      
        // update length
        srparticle.length += ContainedLength(this_point, pos, aa_volumes);

        if (!active_volumes.at(cryostat_index).ContainsPosition(this_point) && active_volumes.at(cryostat_index).ContainsPosition(pos)) {
          exit_point = i-1;
        }

        pos = trajectory.Position(i).Vect();
      }
    }
    if (exit_point < 0 && entry_point >= 0) {
      exit_point = particle.NumberTrajectoryPoints() - 1; 
    }
    if (exit_point >= 0 && ((unsigned)exit_point) < particle.NumberTrajectoryPoints() - 1) {
      srparticle.wallout = GetWallCross(active_volumes.at(cryostat_index), particle.Position(exit_point).Vect(), particle.Position(exit_point+1).Vect()); 
    }

    // other truth information
    srparticle.pdg = particle.PdgCode();

    srparticle.start = (entry_point >= 0) ? particle.Position(entry_point).Vect(): TVector3(-9999, -9999, -9999);
    srparticle.startT = (entry_point >= 0) ? particle.Position(entry_point).T() / 1000. /* ns-> us*/: -9999;
    srparticle.end = (exit_point >= 0) ? particle.Position(exit_point).Vect(): TVector3(-9999, -9999, -9999);
    srparticle.endT = (exit_point >= 0) ? particle.Position(exit_point).T() / 1000. /* ns -> us */ : -9999;
  
    srparticle.startp = (entry_point >= 0) ? particle.Momentum(entry_point).Vect() : TVector3(-9999, -9999, -9999);
    srparticle.startE = (entry_point >= 0) ? particle.Momentum(entry_point).E() : -9999.;
    srparticle.endp = (exit_point >= 0) ? particle.Momentum(exit_point).Vect() : TVector3(-9999, -9999, -9999);
    srparticle.endE = (exit_point >= 0) ? particle.Momentum(exit_point).E() : -9999.;

    srparticle.start_process = GetG4ProcessID(particle.Process());
    srparticle.end_process = GetG4ProcessID(particle.EndProcess());

    srparticle.G4ID = particle.TrackId();

    // See if this MCParticle matches a genie truth
    srparticle.interaction_id = -1;

    art::Ptr<simb::MCTruth> truth = inventory_service.TrackIdToMCTruth_P(particle.TrackId());
    for (unsigned i = 0; i < neutrinos.size(); i++) {
      if (truth.get() == neutrinos[i].get()) {
        srparticle.interaction_id = i;
        break;
      }
    }
  } //FillTrueG4Particle

  void FillFakeReco(const std::vector<art::Ptr<simb::MCTruth>> &mctruths, 
                    const std::vector<art::Ptr<sim::MCTrack>> &mctracks, 
                    const std::vector<geo::BoxBoundedGeo> &volumes,
                    TRandom &rand,
                    std::vector<caf::SRFakeReco> &srfakereco) {
    // iterate and fill
    for (const art::Ptr<simb::MCTruth> mctruth: mctruths) {
      bool do_fill = false;
      caf::SRFakeReco this_fakereco;
      do_fill = FRFillNumuCC(*mctruth, mctracks, volumes, rand, this_fakereco);
      
      // TODO: others?
      // if (!do_fill) ...

      if (do_fill) srfakereco.push_back(this_fakereco);
    }
  }

  std::map<int, std::vector<const sim::IDE*>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels) {
    std::map<int, std::vector<const sim::IDE*>> ret;
    for (const art::Ptr<sim::SimChannel> sc : simchannels) {
      for (const auto &item : sc->TDCIDEMap()) {
        for (const sim::IDE &ide: item.second) {
          // indexing initializes empty vector
          ret[abs(ide.trackID)].push_back(&ide);
        }
      }
    }
    return ret;
  }

} // end namespace


//--------------------------------------------

// helper function definitions

bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance=5.0) {
  TVector3 nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0).Vect();
  TVector3 trkStart = track.Start().Position().Vect();
  return (trkStart - nuVtx).Mag() < distance;
}

// define global static const PDGTable to be used by helper functions
static const TDatabasePDG *PDGTable(new TDatabasePDG);

// returns particle mass in MeV
double PDGMass(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    if (ple == NULL) return -1;
    return ple->Mass() * 1000.0;
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (PDGTable->GetParticle(2212)->Mass() * p +
            PDGTable->GetParticle(2112)->Mass() * n) * 1000.0;
  }
}


bool FRFillNumuCC(const simb::MCTruth &mctruth, 
                  const std::vector<art::Ptr<sim::MCTrack>> &mctracks, 
                  const std::vector<geo::BoxBoundedGeo> &volumes, 
                  TRandom &rand,
                  caf::SRFakeReco &fakereco) {
  // Configuration -- TODO: make configurable?

  // Fiducial Volume
  float xmin = 10.;
  float xmax = 10.;
  float ymin = 10.;
  float ymax = 10.;
  float zmin = 10.;
  float zmax = 100.;

  // energy smearing
  float hadron_smearing = 0.05; 
  float lepton_contained_smearing = 0.02; 
  std::function<float (float)> lepton_exiting_smearing = [](float length) {
    float A = 0.102;
    float B = 0.000612;
    return -A * TMath::Log(B * length);
  };

  // visible energy threshold
  float hadronic_energy_threshold = 0.021; // GeV

  // length cut
  float contained_length_cut = 50.; // cm
  float exiting_length_cut = 100.; // cm

  // first check if neutrino exists
  if (!mctruth.NeutrinoSet()) return false;

  // then check if fiducial
  int cryo_index = -1;
  for (unsigned i = 0; i < volumes.size(); i++) {
    const geo::BoxBoundedGeo &vol = volumes[i];
    geo::BoxBoundedGeo FV(vol.MinX() + xmin, vol.MaxX() - xmax, vol.MinY() + ymin, vol.MaxY() - ymax, vol.MinZ() + zmin, vol.MaxZ() - zmax);
    if (FV.ContainsPosition(mctruth.GetNeutrino().Nu().Position().Vect())) {
      cryo_index = i;
      break;
    }
  }

  if (cryo_index == -1) return false;

  std::vector<geoalgo::AABox> aa_volumes;
  const geo::BoxBoundedGeo &v = volumes.at(cryo_index);
  aa_volumes.emplace_back(v.MinX(), v.MinY(), v.MinZ(), v.MaxX(), v.MaxY(), v.MaxZ());

  // look for CC lepton or a \pi^+/- which can "fake" a numu CC interaction

  int lepton_ind = -1;
  // CC lepton
  if (abs(mctruth.GetNeutrino().Nu().PdgCode()) == 14 && mctruth.GetNeutrino().CCNC() == 0) {
    for (int i = 0; i < (int)mctracks.size(); i++) {
      if (isFromNuVertex(mctruth, *mctracks[i]) && abs(mctracks[i]->PdgCode()) == 13 && mctracks[i]->Process() == "primary") {
        if (lepton_ind == -1 || mctracks[lepton_ind]->Start().E() < mctracks[i]->Start().E()) {
          lepton_ind = i;
        }
      }
    }
  }
  // NC pion
  else if (mctruth.GetNeutrino().CCNC() == 1) {
    for (int i = 0; i < (int)mctracks.size(); i++) {
      if (isFromNuVertex(mctruth, *mctracks[i]) && abs(mctracks[i]->PdgCode()) == 211 && mctracks[i]->Process() == "primary") {
        if (lepton_ind == -1 || mctracks[lepton_ind]->Start().E() < mctracks[i]->Start().E()) {
          lepton_ind = i;
        }
      }
    }
  }

  // no matching track -- no fake reco
  if (lepton_ind == -1) return false;

  // Now set the "lepton"
  caf::SRFakeRecoParticle fake_lepton;
  
  fake_lepton.pid = 13;
  fake_lepton.contained = false;  
  for (const geo::BoxBoundedGeo &vol: volumes) {
    if (vol.ContainsPosition(mctracks[lepton_ind]->Start().Position().Vect()) && vol.ContainsPosition(mctracks[lepton_ind]->Start().Position().Vect())) {
      fake_lepton.contained = true;
    }
  }
  fake_lepton.len = ContainedLength(mctracks[lepton_ind]->Start().Position().Vect(), mctracks[lepton_ind]->End().Position().Vect(), aa_volumes);
  fake_lepton.costh = mctracks[lepton_ind]->Start().Position().Vect().CosTheta();

  // apply length cut
  if (fake_lepton.contained && fake_lepton.len < contained_length_cut) return false;
  if (!fake_lepton.contained && fake_lepton.len < exiting_length_cut) return false;


  // smear the lepton energy
  float smearing = fake_lepton.contained ? lepton_contained_smearing : lepton_exiting_smearing(fake_lepton.len);
  float ke = (mctracks[lepton_ind]->Start().E() - PDGMass(mctracks[lepton_ind]->PdgCode())) / 1000. /* MeV -> GeV*/;
  fake_lepton.ke = rand.Gaus(ke, smearing * ke);
  fake_lepton.ke = std::max(fake_lepton.ke, 0.f);

  // get the hadronic state
  std::vector<caf::SRFakeRecoParticle> hadrons;
  
  for (int i = 0; i < (int)mctracks.size(); i++) {
    if (isFromNuVertex(mctruth, *mctracks[i]) // from this interaction
     && (abs(mctracks[i]->PdgCode()) == 211 || abs(mctracks[i]->PdgCode()) == 321 || abs(mctracks[i]->PdgCode()) == 2212) // hadronic
     && mctracks[i]->Process() == "primary" // primary
     && i != lepton_ind // not the fake lepton
    ) {
      caf::SRFakeRecoParticle hadron; 
      hadron.pid = abs(mctracks[i]->PdgCode());
      hadron.contained = false;
      for (const geo::BoxBoundedGeo &vol: volumes) {
        if (vol.ContainsPosition(mctracks[i]->Start().Position().Vect()) && vol.ContainsPosition(mctracks[i]->Start().Position().Vect())) {
          hadron.contained = true;
        }
      }
      hadron.len = ContainedLength(mctracks[i]->Start().Position().Vect(), mctracks[i]->Start().Position().Vect(), aa_volumes);
      hadron.costh = mctracks[i]->Start().Position().Vect().CosTheta();

      float ke = (mctracks[i]->Start().E() - PDGMass(mctracks[i]->PdgCode())) / 1000. /* MeV -> GeV*/;
      if (ke < hadronic_energy_threshold) continue;

      hadron.ke = rand.Gaus(ke, hadron_smearing * ke);
      hadron.ke = std::max(hadron.ke, 0.f);

      hadrons.push_back(hadron);
    }
  }

  // total up the energy to get the reco neutrino energy
  fakereco.nuE = fake_lepton.ke;   
  for (const caf::SRFakeRecoParticle &had: hadrons) fakereco.nuE += had.ke;
  fakereco.nuE += PDGMass(13) / 1000.; // MeV -.> GeV

  // save the particles
  fakereco.lepton = fake_lepton;
  fakereco.hadrons = hadrons;

  // other info
  TVector3 vertex = mctruth.GetNeutrino().Nu().Position().Vect();
  fakereco.vtx.x = vertex.X(); 
  fakereco.vtx.y = vertex.Y(); 
  fakereco.vtx.z = vertex.Z(); 

  // signal 
  if (abs(mctruth.GetNeutrino().Nu().PdgCode()) == 14 && mctruth.GetNeutrino().CCNC() == 0) {
    fakereco.wgt = 0.8;
  }
  // bkg
  else {
    fakereco.wgt = 1.;
  }
  fakereco.nhad = fakereco.hadrons.size();

  return true;
}

caf::Wall_t GetWallCross(const geo::BoxBoundedGeo &volume, const TVector3 p0, const TVector3 p1) {
  TVector3 direction = (p1 - p0) * ( 1. / (p1 - p0).Mag());
  std::vector<TVector3> intersections = volume.GetIntersections(p0, direction);
  
  assert(intersections.size() == 2);

  // get the intersection point closer to p0
  int intersection_i = ((intersections[0] - p0).Mag() < (intersections[1] - p0).Mag()) ? 0 : 1;

  double eps = 1e-3;
  if (abs(intersections[intersection_i].X() - volume.MinX()) < eps) {
    //std::cout << "Left\n";
    return caf::kWallLeft;
  }
  else if (abs(intersections[intersection_i].X() - volume.MaxX()) < eps) {
    //std::cout << "Right\n";
    return caf::kWallRight;
  }
  else if (abs(intersections[intersection_i].Y() - volume.MinY()) < eps) {
    //std::cout << "Bottom\n";
    return caf::kWallBottom;
  }
  else if (abs(intersections[intersection_i].Y() - volume.MaxY()) < eps) {
    //std::cout << "Top\n";
    return caf::kWallTop;
  }
  else if (abs(intersections[intersection_i].Z() - volume.MinZ()) < eps) {
    //std::cout << "Front\n";
    return caf::kWallFront;
  }
  else if (abs(intersections[intersection_i].Z() - volume.MaxZ()) < eps) {
    //std::cout << "Back\n";
    return caf::kWallBack;
  }
  else assert(false);
  //std::cout << "None\n";  

  return caf::kWallNone;
}//GetWallCross

//------------------------------------------

caf::g4_process_ GetG4ProcessID(const std::string &process_name) {
#define MATCH_PROCESS(name) if (process_name == #name) {return caf::kG4 ## name;}
#define MATCH_PROCESS_NAMED(strname, id) if (process_name == #strname) {return caf::kG4 ## id;}
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
  MATCH_PROCESS_NAMED(xi+Inelastic, xipInelastic)
  MATCH_PROCESS_NAMED(xi-Inelastic, ximInelastic)
  MATCH_PROCESS(kaon0LInelastic)
  MATCH_PROCESS(kaon0SInelastic)
  MATCH_PROCESS(lambdaInelastic)
  MATCH_PROCESS_NAMED(anti-lambdaInelastic, anti_lambdaInelastic)
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
  MATCH_PROCESS(muIoni)
  MATCH_PROCESS(eBrem)
  MATCH_PROCESS(CoulombScat)
  MATCH_PROCESS(hBertiniCaptureAtRest)
  MATCH_PROCESS(hFritiofCaptureAtRest)
  MATCH_PROCESS(photonNuclear)
  MATCH_PROCESS(muonNuclear)
  MATCH_PROCESS(electronNuclear)
  MATCH_PROCESS(positronNuclear)
  MATCH_PROCESS(compt)
  MATCH_PROCESS(eIoni)
  MATCH_PROCESS(muBrems)
  MATCH_PROCESS(hIoni)
  MATCH_PROCESS(muPairProd)
  MATCH_PROCESS(hPairProd)
  MATCH_PROCESS(LArVoxelReadoutScoringProcess)
  std::cerr << "Error: Process name with no match (" << process_name << ")\n";
  assert(false);
  return caf::kG4primary; // unreachable
#undef MATCH_PROCESS
#undef MATCH_PROCESS_NAMED
  
}//GetG4ProcessID
//-------------------------------------------

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

//------------------------------------------------
caf::SRTrackTruth MatchTrack2Truth(const std::vector<art::Ptr<recob::Hit>> &hits) {

  // this id is the same as the mcparticle ID as long as we got it from geant4
  std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(hits, true);
  float total_energy = CAFRecoUtils::TotalHitEnergy(hits);

  caf::SRTrackTruth ret;

  ret.total_deposited_energy = total_energy / 1000. /* MeV -> GeV */;

  // setup the matches
  for (auto const &pair: matches) {
    caf::SRTrackTruth::ParticleMatch match;
    match.G4ID = pair.first;
    match.energy = pair.second / 1000. /* MeV -> GeV */;
    ret.matches.push_back(match);
  }

  // sort highest energy match to lowest
  std::sort(ret.matches.begin(), ret.matches.end(),
	    [](const caf::SRTrackTruth::ParticleMatch &a, const caf::SRTrackTruth::ParticleMatch &b) {
	      return a.energy > b.energy;
	    }
	    );

  return ret;
}//MatchTrack2Truth
//------------------------------------------------
caf::SRTruthMatch MatchSlice2Truth(const std::vector<art::Ptr<recob::Hit>> &hits, 
				   const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
				   const cheat::ParticleInventoryService &inventory_service) {
  caf::SRTruthMatch ret;
  float total_energy = CAFRecoUtils::TotalHitEnergy(hits);
  // speed optimization: if there are no neutrinos, all the matching energy must be cosmic
  if (neutrinos.size() == 0) {
    ret.visEinslc = total_energy / 1000. /* MeV -> GeV */;
    ret.visEcosmic = total_energy / 1000. /* MeV -> GeV */;
    ret.eff = -1; 
    ret.pur = -1;
    ret.index = -1;
    return ret;
  }
  std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(hits, true);
  std::vector<float> matching_energy(neutrinos.size(), 0.);
  for (auto const &pair: matches) {
    art::Ptr<simb::MCTruth> truth;
    try {
      truth = inventory_service.TrackIdToMCTruth_P(pair.first);
    }
    // Ignore track ID's that cannot be looked up
    catch(...) {
      continue;
    }
    for (unsigned ind = 0; ind < neutrinos.size(); ind++) {
      if (truth == neutrinos[ind]) {
	matching_energy[ind] += pair.second;
	break;
      }
    }
  }
  float matching_frac = *std::max_element(matching_energy.begin(), matching_energy.end()) / total_energy;
  int index = (matching_frac > 0.5) ? std::distance(matching_energy.begin(), std::max_element(matching_energy.begin(), matching_energy.end())) : -1;
  float cosmic_energy = total_energy;
  for (float E: matching_energy) cosmic_energy -= E;
  ret.visEinslc = total_energy / 1000. /* MeV -> GeV */;
  ret.visEcosmic = cosmic_energy / 1000. /* MeV -> GeV */;
  ret.index = index;
  if (index >= 0) {
    ret.pur = matching_energy[index] / total_energy;
    // TODO: calculate efficiency 
    ret.eff = 0.;
  }
  else {
    ret.pur = -1;
    ret.eff = -1;
  }
  return ret;
}//Slc2Truth


