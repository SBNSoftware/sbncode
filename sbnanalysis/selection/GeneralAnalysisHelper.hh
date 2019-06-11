#ifndef GENERAL_ANALYSIS_HELPER_H
#define GENERAL_ANALYSIS_HELPER_H

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <utility>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "EventSelectionHelper.hh"

namespace selection{
  
  /**
   * @brief  GeneralAnalysisHelper helper class
   */
  class GeneralAnalysisHelper {

    private : 
      
      /**                                                              
       * @brief  Get the track lengths for a given pdg
       *
       * @param  pdg
       * @param  particle_list
       * @param  lengths
       */
      static void LengthWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &lengths);

      /**                                                              
       * @brief  Get the track opening angles for a given pdg
       *
       * @param  pdg
       * @param  particle_list                                     
       * @param  cos_thetas
       */
      static void CosThetaWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &cos_thetas);
      
      /**                                                              
       * @brief  Returns the energies of particles with a given pdg 
       *
       * @param  pdg
       * @param  particle_list
       * @param  energies
       */
      static void EnergyWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &energies);

      /**                                                              
       * @brief  Returns the Kinetic energies of particles with a given pdg
       *
       * @param  pdg
       * @param  particle_list                                     
       * @param  kinetic_energies
       */
      static void KineticEnergyWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &kinetic_energies);

      /**                                                              
       * @brief  Returns the magnitude of the momentum of particles with a given pdg
       *
       * @param  pdg
       * @param  particle_list                                     
       * @param  momentum_mod
       */
      static void ModulusMomentumWithPdg(const int pdg, const ParticleList &particle_list, std::vector<float> &momentum_mod);
      
      /**
       * @brief  get the MCParticle corresponding to the reco particle general
       *
       * @param  particle reconstructed (track-based) particle
       * @param  mcparticle_list list of mcparticles to loop over
       *
       * @return mcparticle
       *
       */
      static Particle GetMCParticle(const int id, const ParticleList &particle_list );
     
      /**
       * @brief  determine whether the given reconstructed particle's pdg matches the true 
       *
       * @param  particle reconstructed particle to check
       *
       * @return true or false: Matched or not
       *
       */
      static bool MatchedParticle(const Event &e, const Particle &p);
      
      /**
       * @brief  count the number of matched particles with a given pdg code for given particle list
       *
       * @param  event
       * @param  reconstructed particle list
       * @param  pdgcode
       *
       * @return number of matched particles with given pdgcode
       *
       */
      static unsigned int CountMatchedParticles(const Event &e, const ParticleList &particle_list, const int pdg);

    public : 

      typedef std::vector<Particle> ParticleList;
      typedef std::vector<Event>    EventList;

      /**
       * @brief  Get NuMu topology map
       */
      static TopologyMap GetNuMuTopologyMap();

      /**
       * @brief  Get NC topology map
       */
      static TopologyMap GetNCTopologyMap();

      /**
       * @brief  Get NC 0Pi topology map
       */
      static TopologyMap GetNC0PiTopologyMap();

      /**
       * @brief  Get NC 1Pi topology map
       */
      static TopologyMap GetNC1PiTopologyMap();

      /**
       * @brief  Get NC 2Pi topology map
       */
      static TopologyMap GetNC2PiTopologyMap();

      /**
       * @brief  Get CC inclusive topology map
       */
      static TopologyMap GetCCIncTopologyMap();

      /**
       * @brief  Get CC 0Pi topology map
       */
      static TopologyMap GetCC0PiTopologyMap();

      /**
       * @brief  Get CC 0Pi 1Proton topology map
       */
      static TopologyMap GetCC0Pi1PTopologyMap();

      /**
       * @brief  Get CC 0Pi 2Protons topology map
       */
      static TopologyMap GetCC0Pi2PTopologyMap();

      /**
       * @brief  Get CC 0Pi 3Protons topology map
       */
      static TopologyMap GetCC0Pi3PTopologyMap();

      /**
       * @brief  Get CC 0Pi 5Protons topology map
       */
      static TopologyMap GetCC0Pi5PTopologyMap();

      /**
       * @brief  Get CC 1Pi topology map
       */
      static TopologyMap GetCC1PiTopologyMap();

      /**
       * @brief  Get CC 2Pi topology map
       */
      static TopologyMap GetCC2PiTopologyMap();

      /**
       * @brief  Get CC 1Pi0 topology map
       */
      static TopologyMap GetCCPi0TopologyMap();

      /*8
       * @brief  Get the NuE topology map
       */
      static TopologyMap GetNuETopologyMap();
      /**
       * @brief  Get the number of escaping reconstructed tracks or MCParticles
       *
       * @param  e Current event
       *
       * @return Number of escaping reconstructed tracks or MCParticles
       */
      static unsigned int NumberEscapingTracks(const Event &e);
      
      /**
       * @brief  Finds if there is more than one escaping track in an event
       *
       * @param  e Current event
       *
       * @return True if there is a maximum of one, false if there are more than one
       */
      static bool MaxOneEscapingTrack(const Event &e);
      
      /**                                                              
       * @brief  Gives the number of MC, Reco and Coincidences for a given topology                                                           
       *
       * @param  e current event
       * @param  signal_map_topology chosen topology
       * @param  count_true number of true events with chosen topology
       * @param  count_signal number of signal events with chosen topology
       * @param  count selected number of selected events with chosen topology
       */
      static void TopologyStatistics(const Event &e, const TopologyMap signal_map_topology, double &count_true, double &count_signal, double &count_selected);
     
      /**                                                              
       * @brief  Obtains the Topology matrix for a specific set of events                                                                    
       *
       * @param  e current event
       * @param  count_true_topology number of true events with chosen topology
       * @param  count_signal_topology number of signal events with chosen topology
       * @param  count selected_topology number of selected events with chosen topology
       *
       * @return Matrix of topology-based statistics
       */
      static ParticleMatrix TopologyMatrix(const Event &e, ParticleMatrix &count_true_topology, ParticleMatrix &count_signal_topology, ParticleMatrix &count_selected_topology);
  
      /**                                                              
       * @brief  Get the longest MC tracks with a given pdg     
       *
       * @param  event
       * @param  pdg                                                 
       * @param  lengths 
       */
      static void GetMCLengthWithPdg(const Event &e, const int pdg, std::vector<float> &lengths);

      /**                                                              
       * @brief  Get the longest reconstructed tracks with a given pdg
       *
       * @param  event
       * @param  pdg                                                    
       * @param  lengths 
       */
      static void GetRecoLengthWithPdg(const Event &e, const int pdg, std::vector<float> &lengths);

      /**                                                              
       * @brief  Get the cos thetas for given MC pdg
       *
       * @param  event
       * @param  pdg
       * @param  cos_thetas
       */
      static void GetMCCosThetaWithPdg(const Event &e, const int pdg, std::vector<float> &cos_thetas);

      /**                                                              
       * @brief  Get the cos theta with reco pdg
       *
       * @param  event
       * @param  pdg
       * @param  cos_thetas
       */
      static void GetRecoCosThetaWithPdg(const Event &e, const int pdg, std::vector<float> &cos_thetas);
      
      /**                                                              
       * @brief  Get the MC energies with a given pdg
       *
       * @param  event
       * @param  pdg
       * @param  energies
       */
      static void GetMCEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &energies);

      /**                                                              
       * @brief  Get the reco energy with pdg
       *
       * @param  event
       * @param  pdg
       * @param  energies
       */
      static void GetRecoEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &energies);

      /**                                                              
       * @brief  Get the MC kinetic energy with pdg
       *
       * @param  event
       * @param  pdg
       * @param  kinetic_energies
       */
      static void GetMCKineticEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &kinetic_energies);

      /**                                                              
       * @brief  Get the reco kinetic energy with pdg
       *
       * @param  event
       * @param  pdg
       * @param  kinetic_energies
       */
      static void GetRecoKineticEnergyWithPdg(const Event &e, const int pdg, std::vector<float> &kinetic_energies);
      
      /**                                                              
       * @brief  Get the MC modulus of the momentum with pdg
       *
       * @param  event
       * @param  pdg
       * @param  momentum_mod
       */
      static void GetMCModulusMomentumWithPdg(const Event &e, const int pdg, std::vector<float> &momentum_mod);
      
      /**                                                              
       * @brief  Get the reco modulus of the momentum with pdg
       *
       * @param  event
       * @param  pdg
       * @param  momentum_mod
       */
      static void GetRecoModulusMomentumWithPdg(const Event &e, const int pdg, std::vector<float> &momentum_mod);
      
      /**
       * @brief  Calculates the Efficiency, Purity, Background Rejection
       *         Parameters for Efficiency calculation ( MC, signal and selected ) for a given topology : 
       *         0-> No muon, 
       *         1 -> CCinclusive,
       *         2-> CC0pi, 
       *         3-> CC1pi+/-,
       *         4-> CC1pi0
       *
       * @param  count_mc
       * @param  count_signal 
       * @param  count_selected 
       **/
      static double Efficiency(const std::vector< double > &count_mc, const std::vector< double > &count_signal, const std::vector< double > &count_selected, const TopologyMap &topology);

      /**
       * @brief  Save Topology Matrix into a file
       *         BackGround Study : topology mis identification table 
       *         0-> No muon, 
       *         1 -> CCinclusive,
       *         2-> CC0pi, 
       *         3-> CC1pi+/-,
       *         4-> CC1pi0
       *
       * @param  count_mc_topology 
       * @param  count_signal_topology 
       * @param  count_selected_topology 
       **/
      static void SaveTopologyMatrix(const ParticleMatrix &count_mc_topology, const ParticleMatrix &count_signal_topology, const ParticleMatrix &count_selected_topology);

      /**
       * @brief Saves the event information in a file ( types of particles in the event 
       * and topology for True and Selected
       **/
      static void EventInformationParticles(const Event &e, const std::string name, const int event_number);

      /**
       * @brief Saves the event characteristics in a file ( length , angle and kinetic energy ) 
       * for the selected topology
       **/
      static void EventProperties(const Event &e, const TopologyMap &topology, std::string event_file, const int event_number);
      
      /**
       * @brief  get the MCParticle corresponding to the reco particle via charge
       *
       * @param  event
       * @param  particle reconstructed (track-based) particle
       *
       * @return mcparticle
       *
       */
      static Particle GetMCParticleCharge(const Event &e, const Particle &particle);

      /**
       * @brief  get the MCParticle corresponding to the reco particle via energy
       *
       * @param  event
       * @param  particle reconstructed (track-based) particle
       *
       * @return mcparticle
       *
       */
      static Particle GetMCParticleEnergy(const Event &e, const Particle &particle);

      /**
       * @brief  get the MCParticle corresponding to the reco particle via hits
       *
       * @param  event
       * @param  particle reconstructed (track-based) particle
       *
       * @return mcparticle
       *
       */
      static Particle GetMCParticleHits(const Event &e, const Particle &particle);
  
      /**
       * @brief  If the reconstructed particle has a match by
       *            0 == hits
       *            1 == charge
       *            2 == energy
       *          get the MCParticle using the relevant method in order of preference (0 -> 2)
       *
       * @param  event
       * @param  particle reconstructed (track-based) particle
       *
       * @return mcparticle
       *
       */
      static Particle GetBestMCParticle(const Event &e, const Particle &particle);
      
      /**
       * @brief  from a selected event with a given topology, count the number of MC particles
       * with a given pdg code
       *
       * @param  event
       * @param  topology
       * @param  pdgcode
       *
       * @return number of MC particles for selected event of a given pdgcode
       */
      static unsigned int CountMCParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg);

      /**
       * @brief  from a signal event with a given topology, count the number of MC particles
       * with a given pdg code
       *
       * @param  event
       * @param  topology
       * @param  pdgcode
       *
       * @return number of MC particles for signal event of a given pdgcode
       */
      static unsigned int CountMCParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg);

      /**
       * @brief  from a selected event with a given topology, count the number of Reco particles
       * with a given pdg code
       *
       * @param  event
       * @param  topology
       * @param  pdgcode
       *
       * @return number of Reco particles for selected event of a given pdgcode
       */
      static unsigned int CountRecoParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg);

      /**
       * @brief  from a signal event with a given topology, count the number of Reco particles
       * with a given pdg code
       *
       * @param  event
       * @param  topology
       * @param  pdgcode
       *
       * @return number of Reco particles for signal event of a given pdgcode
       */
      static unsigned int CountRecoParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg);

      /**
       * @brief  from a selected event with a given topology, count the number of matched particles
       * with a given pdg code
       *
       * @param  event
       * @param  topology
       * @param  pdgcode
       *
       * @return number of matched particles for selected event of a given pdgcode
       */
      static unsigned int CountMatchedParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int pdg);

      /**
       * @brief  from a signal event with a given topology, count the number of matched particles
       * with a given pdg code
       *
       * @param  event
       * @param  topology
       * @param  pdgcode
       *
       * @return number of matched particles for signal event of a given pdgcode
       */
      static unsigned int CountMatchedParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int pdg);
      
      /**
       * @brief  count the number of matched particles with a given pdg code for whole event
       *
       * @param  event
       * @param  pdgcode
       *
       * @return number of matched particles with given pdgcode
       */
      static unsigned int CountMatchedParticlesAll(const Event &e, const int pdg);

      /**
       * @brief  Count mismatched particles based on topology and given pdg codes
       *
       * @param  event
       * @param  topology
       * @param  true pdg code of particle
       * @param  reconstructed pdg code of particle
       *
       * @return number of mis-matched particles for given pdg codes selected events
       */
      static unsigned int CountMisMatchedParticlesByTopologySelected(const Event &e, const TopologyMap &topology, const int true_pdg, const int reco_pdg);
      
      /**
       * @brief  Count mismatched particles based on topology and given pdg codes
       *
       * @param  event
       * @param  topology
       * @param  true pdg code of particle
       * @param  reconstructed pdg code of particle
       *
       * @return number of mis-matched particles for given pdg codes
       */
      static unsigned int CountMisMatchedParticlesByTopologySignal(const Event &e, const TopologyMap &topology, const int true_pdg, const int reco_pdg);
      
      /**
       * @brief  Count mismatched particles and given pdg codes
       *
       * @param  event
       * @param  true pdg code of particle
       * @param  reconstructed pdg code of particle
       *
       * @return number of mis-matched particles for given pdg codes
       */
      static unsigned int CountMisMatchedParticles(const Event &e, const int true_pdg, const int reco_pdg);
      
      /**
       * @brief  for all events, count all matched particles by topology and fill file
       *
       * @param  event list
       * @param  topology
       * @param  topology name
       * @param  file to append
       *
       */
      static void FillTopologyBasedParticleStatisticsFile(const EventList &ev_list, const TopologyMap &topology, const std::string &topology_name, std::ofstream &os);

      /**
       * @brief  for all events, count all mismatched particles by topology and fill file
       *
       * @param  event list
       * @param  topology
       * @param  topology name
       * @param  file to append
       *
       */
      static void FillTopologyBasedParticleMisIdStatisticsFile(const EventList &ev_list, const TopologyMap &topology, const std::string &topology_name, std::ofstream &os);

      /**
       * @brief  for all events, count all matched particles and fill file
       *
       * @param  event list
       * @param  file name
       *
       */
      static void FillGeneralParticleStatisticsFile(const EventList &ev_list, std::ofstream &os);

      /**
       * @brief  for all events, count all mismatched particles and fill file
       *
       * @param  event list
       * @param  file name
       *
       */
      static void FillGeneralParticleMisIdStatisticsFile(const EventList &ev_list, std::ofstream &os);

      /**
       * @brief  find out if a reconstructed particle has a matching truth particle
       *          first check using hits, return 0 if true
       *          second check using charge, return 1 if true
       *          third check using energy, return 2 if true
       *          if no match, return -1
       *
       * @param  current event
       * @param  reconstructed particle
       *
       * @return enumerated result
       *
       */
      static int ParticleHasAMatch(const Event &e, const Particle &p);

      /**
       * @brief  Find out if an MC particle has a corresponding reconstructed particle
       *
       * @param  event
       * @param  mc particle
       *
       * @return true or false: true == has reconstructed particle, false == hasn't
       *
       */
      static bool HasBeenReconstructed(const Event &e, const Particle &p);
      

  }; // GeneralAnalysisHelper
} // namespace: selection
#endif
