#ifndef SBNCODE_SINGLEPHOTONANALYSIS_RECO_TRUTH_MATCHING_H
#define SBNCODE_SINGLEPHOTONANALYSIS_RECO_TRUTH_MATCHING_H

#include <vector>
#include <map>

#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"


//#include <typeinfo>

namespace single_photon
{
  void RecoMCTracks(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Track>>& tracks,  
      std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,
      std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,
      std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector,  
      std::map< int, art::Ptr<simb::MCParticle> > &      MCParticleToTrackIdMap, 
      std::vector<double> & vfrac
      );

  //recoMCmatching but specifically for recob::showers
  void showerRecoMCmatching(
      std::vector<PandoraPFParticle> all_PPFPs,
      std::vector<art::Ptr<recob::Shower>>& showerVector,
      std::map<art::Ptr<recob::Shower>,art::Ptr<simb::MCParticle>>& showerToMCParticleMap,
      art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
      std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
      std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap);



  /* @brief: a simpler MCmatching function for track and shower
   * @argument to be filled in function body:
   *     objectToMCParticleMap: map of object (track, shower) to its best-matching MCParticle
   *     mcParticleVector: a vector of best-matching MCParticle corresponding to objectVector
   * @return: a vector of fraction number, which is the fraction of unassociated hits in all reco hits of PFParticle
   */  
  //Typenamed for recob::Track and recob::Shower
    std::vector<double> trackRecoMCmatching(std::vector<art::Ptr<recob::Track>>& objectVector,
        std::map<art::Ptr<recob::Track>,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
        std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
        art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
        std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector);


  int    photoNuclearTesting(std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector);

}//namespace end

#endif // SBNCODE_SINGLEPHOTONANALYSIS_RECO_TRUTH_MATCHING_H
