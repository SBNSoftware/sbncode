#ifndef HELPER_PANDORAPFPARTICLES_H
#define HELPER_PANDORAPFPARTICLES_H


#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "sbnobj/Common/CRT/CRTHit.hh"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "messagefacility/MessageLogger/MessageLogger.h"


#include <utility>
#include <iostream>
#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <sys/stat.h>



namespace single_photon
{
  class PandoraPFParticle{

    private:
    double pVertex_pos[3] = {-9999,-9999,-9999};

    double pNuScore = -999 ;
    double pTrackScore = -999;

    bool pIsNeutrino = false;
    bool pIsClearCosmic = false;
    bool pIsNuSlice = false; //it is defined in the DefineNuSlice() function.
    bool pHasPID = false;

    int pHasShower = 0;
    int pHasTrack = 0;
    int pPdgCode = -999;
    int pPFParticleID = -9;
    int pAncestorID = -9;
    int pSliceID = -9;


  public:
    //constructor:
    PandoraPFParticle( 
        art::Ptr<recob::PFParticle> input_PFParticle,
        std::vector< art::Ptr< larpandoraobj::PFParticleMetadata > > input_MetaData,
        std::vector< art::Ptr<recob::Vertex > > input_Vertex,
        std::vector< art::Ptr<recob::Cluster> > input_Clusters,
        std::vector< art::Ptr<recob::Shower > > input_Showers,
        std::vector< art::Ptr<recob::Track  > > input_Tracks,
        art::FindManyP<recob::Hit>   input_Hits );

    art::Ptr< recob::PFParticle > pPFParticle;
    art::Ptr< recob::Shower>  pShower; //with 0 or 1 element
    art::Ptr< recob::Track >  pTrack; //with 0 or 1 element
    art::Ptr< recob::Slice >  pSlice; //in helper_connector.h get the id from pSlice->key()
    art::Ptr< recob::PFParticle > pAncestor; //found by tracing Parent()

    art::Ptr<anab::ParticleID>  pParticleID; //for track only;
    art::Ptr< simb::MCTruth >  pMCTruth;//WARNING NOT USED YET

    std::vector< art::Ptr< larpandoraobj::PFParticleMetadata > > pMetaData;
    std::vector< art::Ptr< recob::Vertex > > pVertex;
    std::vector< art::Ptr< recob::Hit > >  pSliceHits;
    std::vector< art::Ptr< recob::Hit > >  pPFPHits;
    std::vector< art::Ptr< recob::Cluster > > pClusters;
    std::vector<art::Ptr<anab::Calorimetry>> pCalorimetries;
    std::vector< art::Ptr< recob::SpacePoint > > pSpacePoints;//WARNING NOT USED YET
    std::vector< art::Ptr< simb::MCParticle > > pMCParticles;//WARNING NOT USED YET


    //set methods
    void set_NuScore (const double input_score){ pNuScore = input_score; }
    void set_IsNuSlice (const bool input_bool){ pIsNuSlice = input_bool; }
    void set_HasPID (const bool input_bool){ pHasPID = input_bool; }

    void set_AncestorID (const int input_number){ pAncestorID = input_number; }
    void set_SliceID (const int input_number){ pSliceID = input_number; }

    void set_ParticleID (const art::Ptr<anab::ParticleID>  input_ParticleID){pParticleID = input_ParticleID; }
    void set_Calorimetries (const std::vector<art::Ptr<anab::Calorimetry>>  input_Calorimetries){pCalorimetries = input_Calorimetries; }

    //call method
    const art::Ptr<anab::ParticleID>  get_ParticleID() const;
    const std::vector<art::Ptr<anab::Calorimetry>> get_Calorimetries() const;


    const double* get_Vertex_pos()     const;

    const double  get_NuScore()        const;
    const double  get_TrackScore ()    const;

    const bool    get_IsNeutrino ()    const;
    const bool    get_IsClearCosmic () const;
    const bool    get_IsNuSlice ()     const;
    const bool    get_HasPID ()        const;

    const int     get_HasShower ()     const;
    const int     get_HasTrack ()      const;
    const int     get_PdgCode()        const;
    const int     get_PFParticleID()   const;
    const int     get_AncestorID()     const;
    const int     get_SliceID ()       const;
  };

  inline const art::Ptr<anab::ParticleID> PandoraPFParticle::get_ParticleID() const { return pParticleID; }
  inline const std::vector<art::Ptr<anab::Calorimetry>> PandoraPFParticle::get_Calorimetries() const { return pCalorimetries; }

  inline const double* PandoraPFParticle::get_Vertex_pos()    const { return pVertex_pos; }

  inline const double PandoraPFParticle::get_NuScore()        const { return pNuScore; }
  inline const double PandoraPFParticle::get_TrackScore ()    const { return pTrackScore; }

  inline const bool   PandoraPFParticle::get_IsNeutrino ()    const { return pIsNeutrino; }
  inline const bool   PandoraPFParticle::get_IsClearCosmic () const { return pIsClearCosmic; }
  inline const bool   PandoraPFParticle::get_IsNuSlice ()     const { return pIsNuSlice; }
  inline const bool   PandoraPFParticle::get_HasPID ()     const { return pHasPID; }

  inline const int    PandoraPFParticle::get_HasShower ()     const { return pHasShower; }
  inline const int    PandoraPFParticle::get_HasTrack ()      const { return pHasTrack; }
  inline const int    PandoraPFParticle::get_PdgCode()        const { return pPdgCode; }
  inline const int    PandoraPFParticle::get_PFParticleID()   const { return pPFParticleID; }
  inline const int    PandoraPFParticle::get_AncestorID()     const { return pAncestorID; }
  inline const int    PandoraPFParticle::get_SliceID ()       const { return pSliceID; }



//helper functions exclusively for PandoraPFParticles

//1. Filler for each PandoraPFParticle
  //find pAncestor & pAncesotrID
  void PPFP_FindAncestor ( std::vector< PandoraPFParticle > & PPFPs);

  //Fill pSlice, pHits, pSliceID
  void PPFP_FindSliceIDandHits(std::vector< PandoraPFParticle > & PPFPs, art::Ptr< recob::Slice >  slice, const std::vector<art::Ptr<recob::PFParticle> > PFP_in_slice, const std::vector<art::Ptr<recob::Hit> > Hit_inslice);

  //refill pNuScore and pIsNuSlice
  int DefineNuSlice(std::vector< PandoraPFParticle > & PPFPs);

//2. Trackers to find the correct PandoraPFParticle based on different inputs
  PandoraPFParticle *PPFP_GetPPFPFromShower( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Shower> pShower);
  PandoraPFParticle* PPFP_GetPPFPFromTrack( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Track> pTrack);
  PandoraPFParticle* PPFP_GetPPFPFromPFID( std::vector< PandoraPFParticle > & PPFPs, int id);

}

#endif
