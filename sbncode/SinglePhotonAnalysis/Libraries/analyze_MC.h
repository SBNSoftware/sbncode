#ifndef SBNCODE_SINGLEPHOTONANALYSIS_ANALYZE_MC_H
#define SBNCODE_SINGLEPHOTONANALYSIS_ANALYZE_MC_H

#include <vector>
#include <map>

#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"

namespace single_photon
{

//analyze_Geant4.h
  void AnalyzeGeant4( const    std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector);

//analyze_EventWeight.h
  void AnalyzeEventWeight(art::Event const & e);

//analyze_MCTruth.h
  void AnalyzeRecoMCSlices(std::string signal_def, 
      std::vector<PandoraPFParticle> all_PPFPs,
      std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,
      std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
      std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap);

  //This only look at MCTruch info. Reco matching create sim_shower/track for pairing up MCTruth to Reco objects;
  void AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector , std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector);


}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_ANALYZE_MC_H
