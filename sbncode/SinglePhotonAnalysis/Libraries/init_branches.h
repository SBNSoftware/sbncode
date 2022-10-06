#ifndef SBNCODE_SINGLEPHOTONANALYS_INIT_BRANCHES_H
#define SBNCODE_SINGLEPHOTONANALYS_INIT_BRANCHES_H

#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"

#include "art/Framework/Principal/Event.h"

#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"

namespace single_photon
{

  bool g_is_verbose = true;

  /**
   * @brief: reset/clear data members
   */
  void ClearMeta(var_all& vars);
  void CreateMetaBranches(var_all& vars);

  void ClearStubs(var_all& vars);
  void CreateStubBranches(var_all& vars);

  void ClearSecondShowers(var_all& vars); /* reset and clear variable/vectors related to second shower */
//  void ResizeSecondShowers(size_t size); /* currently does nothing */ 
  void CreateSecondShowerBranches(var_all& vars); /*create branches in vertex tree for second shower related variables */

  void ClearSecondShowers3D(var_all& vars); /* reset and clear variables/vectors realted to second shower 3D */
  void CreateSecondShowerBranches3D(var_all& vars); /*create branches in vertex tree for second shower 3D */


  void ClearIsolation(var_all& vars);  /* clear vector members related to isolation */
  void CreateIsolationBranches(var_all& vars);  /* create branches for vectors related to isolation in vertex_tree */

  void ClearFlashes(var_all& vars);  /* clear and reset all the flash-related vectors/variables */
  void ResizeFlashes(size_t size, var_all& vars); /* resize flash-related vectors */
  void CreateFlashBranches(var_all& vars); /* create branches for flashes in vertex_tree */

  void ClearTracks(var_all& vars);  /* clear track related variable and vectors */
  void ResizeTracks(size_t size, var_all& vars);  /* resize track related vectors */
  void CreateTrackBranches(var_all& vars); /* create track related branch in vertex tree */

  void ClearShowers(var_all& vars);  /* clear and reset shower related vectors/variables, for reco, and sim shower */
  void ResizeShowers(size_t size, var_all& vars);   /* resize vectors related to reco, and sim showers */
  void CreateShowerBranches(var_all& vars);  /* create branches for shower-related vec/vars in vertex_tree */

  void ClearMCTruths(var_all& vars);
  void ResizeMCTruths(size_t size, var_all& vars);  /* resize mctruth daughters vectors */
  void CreateMCTruthBranches(var_all& vars);

  /**
   * @brief: fill event weight related variables */
  void AnalyzeEventWeight(art::Event const & e );
  void ClearEventWeightBranches(var_all& vars);  /* reset eventweight related variable */
  void CreateEventWeightBranches(var_all& vars);  /* create branches for eventweight related variable in eventweight_tree */


  /**
   * @brief: fill event weight related variables */
  void ClearGeant4Branches(var_all& vars);  /* reset eventweight related variable */
  void CreateGeant4Branches(var_all& vars);  /* create branches for eventweight related variable in eventweight_tree */
  void AnalyzeGeant4( const    std::vector<art::Ptr<simb::MCParticle>> &mcParticleVector,var_all& vars);    


  void ClearSlices(var_all& vars); /* reset and clear variables/vectors related to slice */
  void ResizeSlices(size_t size, var_all& vars);   /* resize vectors related to slice */ 
  void CreateSliceBranches(var_all& vars);  /* create slice branches in ncdelta_slice_tree and vertex_tree */


  void Save_EventMeta( art::Event &evt, var_all& vars);
  void Save_PFParticleInfo( std::vector<PandoraPFParticle> PPFPs, var_all& vars, para_all& paras);

}
#endif // SBNCODE_SINGLEPHOTONANALYS_INIT_BRANCHES_H
