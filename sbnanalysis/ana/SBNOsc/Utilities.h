#ifndef __sbnanalysis_ana_SBNOsc_Utilities__
#define __sbnanalysis_ana_SBNOsc_Utilities__

/**
 * \file Utilities.h
 *
 * Common utilties
 *
 * This is some auxiliary code that is not a selection, but does a piece
 * of the analysis. We can define any number of other functions, classes,
 * etc. which we use in the selection.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>
 */

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "gallery/Event.h"

#include "core/Event.hh"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"

#include "TRandom3.h"

namespace ana {
  namespace SBNOsc {

/** A function that says hello. */
void hello();


/** Extract truth information to approximate reconstruction. */
event::Interaction TruthReco(const simb::MCTruth& mctruth);


/**
 * Get oscillation probability of muon neutrino in a 3+1 model. I.e.
 * probability that the numu will stay a numu.
 *
 * \param numu_energy Energy of incident muon neutrino in GeV
 * \param numu_dist Distance travelled by muon neutrino in km
 * \param osc_dm2 dm^2 of sterile netrino in eV^2
 * \param osc_angle Sterile neutrino mixing angle
 *
 * \return Probability of muon neutrino not oscillating in 3+1 model.
 */
double NuMuOscillation(double numu_energy, double numu_dis,
                       double osc_dm2, double osc_angle);


/**
 * Finds length of line segment contained inside AABox. Make sure that
 * AABox and TVector's use the same units.
 *
 * \param v0 the first point of the line segment
 * \param v1 the second point of the line segment
 * \param boxes a list of fiducial volumes instantiated as AABoxes
 *
 * \return Length of line segment contained in the list of AABox's.
 */
double containedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes);


/**
 * Get mass from PDGID of particle in MeV/c^2.
 *
 * \param pdg The Particle Data Group ID of the particle (as returned
 *        by i.e. an MCTruth object)
 *
 * \return Mass of particle in MeV/c^2
 */
double PDGMass(int pdg);


/**
 * Get charge from PDGID of particle in |e|/3.
 * \param pdg The Particle Data Group ID of the particle (as returned by
 *        i.e. an MCTruth object)
 *
 * \return Charge of particle in |e|/3.
 */
double PDGCharge(int pdg);


/**
 * Returns whether track/shower object is from the neutrino vertex
 *
 * \param mc MCTruth corresponding to neutrino interaction
 * \param show The object to be matched
 * \param distance between shower start and interaction vertex
 * \return Whether track/shower object is from neutrino vertex
 */
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                           float distance=5.0);


/**
 * Returns whether track/shower object is from the neutrino vertex
 *
 * \param mc MCTruth corresponding to neutrino interaction
 * \param track The object to be matched
 * \param distance between track start and interaction vertex
 * \return Whether track/shower object is from neutrino vertex
 */
bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                            float distance=5.0);

/**
 * Returns whether a mc particle object is from a netrino vertex
 * \param mc MCTruth corresponding to neutrino interaction
 * \param particle The object to be matched
 * \param distance between track start and interaction vertex
 * \return Whether track/shower object is from neutrino vertex
 */
  
bool isFromNuVertex(const simb::MCTruth& mc, const simb::MCParticle* &particle, float distance=5.0);

 bool isFromNuVertex(std::vector<simb::MCTruth>& mcs,std::map<int, const simb::MCParticle*>& mcparticles, int particle_id, float distance=5.0);

bool isChargedPrimary(const simb::MCTruth& mc, std::map<int, const simb::MCParticle*>& mcparticles, int particle_id);

/**
 * Calculate CCQE energy from associated lepton information (and optional
 * distortion). Energy in GeV.
 *
 * \param l_momentum Lepton momentum (in any units -- used only to
 *        get angle info)
 * \param l_energy Lepton energy in GeV
 * \param energy_distortion Optional energy distortion in GeV
 * \param angle_distortion Optiona langle distortion
 *
 * \return CCQE energy in GeV.
 */
double ECCQE(const TVector3& l_momentum, double l_energy,
             double energy_distortion=0., double angle_distortion=0.);


/** 
 * \class Struct containing information used in calculation of visible Energy
 *
 */
struct VisibleEnergyCalculator {
  int lepton_pdgid; //!< PDGID of lepton in the interaction. Used to add 
    // in energy corresponding to lepton mass for CC events (and confused NC
    // events). If you don't want to add in a lepton mass to the energy
    // accounting, set it to 0. 
  double track_threshold; //!< Energy threshold of track energy counted in calculation [GeV].
  double shower_threshold; //!< Energy threshold of shower energy counted in calculation [GeV].
  double track_energy_distortion; //!< Distortion of energies of tracks (%).
  double shower_energy_distortion; //!< Distortion of energies of showers (%).
    
  double lepton_energy_distortion_contained; //!< Distortion of energies of primary lepton whose tracks are contained within the TPC (%).
  double lepton_energy_distortion_leaving_A; //!< Parameter in function to calculate primary lepton energy resolution. 
      // (%) = -A * Log(B * L)  where L is the lepton contained length
  double lepton_energy_distortion_leaving_B; //!< Parameter in function to calculate primary lepton energy resolution.
      // (%) = -A * Log(B * L)  where L is the lepton contained length
  int lepton_index; //!< Index of lepton in the mctrack object
  bool lepton_contained; //!< True if primary lepton's track is contained within TPC.
  double lepton_contained_length; //!< Length of section of primary lepton's track that is contained within the TPC.

  VisibleEnergyCalculator(): 
    lepton_pdgid(0),
    track_threshold(0),
    shower_threshold(0),
    track_energy_distortion(0),
    shower_energy_distortion(0),
    lepton_energy_distortion_contained(0),
    lepton_energy_distortion_leaving_A(0),
    lepton_energy_distortion_leaving_B(0),
    lepton_contained(false),
    lepton_contained_length(1)
  {}
};


/**
 * Get the "visible" energy from a neutrino interaction. Is equal to sum of
 * non-neutral hadronic kinetic energies and lepton total energies.
 *
 * \param ev The gallery event.
 * \param mctruth The MCTruth object corresponding to the interaction.
 * \param mctrack_list Vector of MCTrack objects in the gallery event.
 * \param mcshower_list Vector of MCShower objects in the gallery event.
 * \param calculator Struct containing values to be used in energy calculation
 * \param smeared_lepton_energy lepton energy to be used in calculation -- will default to smearLeptonEnergy(mctruth, calculator) if not set
 *
 * \return Visble energy in GeV.
 * */
 double visibleEnergy(TRandom3& rand, const simb::MCTruth &mctruth, const std::vector<sim::MCTrack> &mctrack_list, const std::vector<sim::MCShower> &mcshower_list,  
		     const VisibleEnergyCalculator &calculator=VisibleEnergyCalculator(), bool include_showers=true);

/** 
 * Get the seperate hadronic and leptonic energy from the event. The first element in the vector is the hadronic energy the second is the leptonic. 
 * */
 std::vector<double> FlavourEnergyDeposition(TRandom3& rand, const simb::MCTruth &mctruth, std::map<int,const simb::MCParticle*>& mcparticles,std::map<int,double>& mcvisibleparticles, std::vector<geoalgo::AABox>& Volumes, const VisibleEnergyCalculator &calculator=VisibleEnergyCalculator());

/** Get the smeared energy from a lepton.
 * \param mctrack The MCTrack object corresponding to the lepton
 * \param calculator Struct containing values to be used in energy calculation
 *
 * */
 double smearLeptonEnergy(TRandom3& rand,const sim::MCTrack &mct, const VisibleEnergyCalculator &calculator=VisibleEnergyCalculator());

/** Get the smeared energy from a lepton.
 * \param lepton The MCParticle object object corresponding to the lepton
 * \param calculator Struct containing values to be used in energy calculation
 *
 * */
 double smearLeptonEnergy(TRandom3& rand, const simb::MCParticle* &lepton, const VisibleEnergyCalculator &calculator=VisibleEnergyCalculator());

 double smearLeptonEnergy(TRandom3& rand, double& LeptonE, int& PdgCode,  const VisibleEnergyCalculator &calculator=VisibleEnergyCalculator());

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_Utilities__
