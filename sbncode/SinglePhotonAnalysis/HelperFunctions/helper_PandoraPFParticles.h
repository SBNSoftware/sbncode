#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

//KENG, should use these in sbndcode

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
//KENG
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include  "nusimdata/SimulationBase/GTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

//#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//#include "ubobj/Optical/UbooneOpticalFilter.h"

// Helper function for PID stuff
//#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TRandom3.h"
#include "TGeoPolygon.h"

//#include "Pandora/PdgTable.h"
#include <chrono>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <vector>
#include <set>
#include <sys/stat.h>



namespace single_photon
{
  class PandoraPFParticle{

    private:
      double pVertex_pos[3] = {-9999,-9999,-9999};//d

    double pNuScore = -999 ;//d
    double pTrackScore = -999;//d

    bool pIsNeutrino = false;//d
    bool pIsClearCosmic = false;//d
    bool pIsNuSlice = false;//d it is defined in the DefineNuSlice() function.
    bool pHasPID = false;//d

    int pHasShower = 0;//d
    int pHasTrack = 0;//d
    int pPdgCode = -999;//d
    int pPFParticleID = -9;//d
    int pAncestorID = -9;//d
    int pSliceID = -9;//d


  public:
    //constructor:
    PandoraPFParticle( 
        art::Ptr<recob::PFParticle> input_PFParticle,
        std::vector< art::Ptr< larpandoraobj::PFParticleMetadata > > input_MetaData,
        std::vector< art::Ptr<recob::Vertex > > input_Vertex,
        std::vector< art::Ptr<recob::Cluster> > input_Clusters,
        std::vector< art::Ptr<recob::Shower > > input_Showers,
        std::vector< art::Ptr<recob::Track  > > input_Tracks,
        art::FindManyP<recob::Hit>   input_Hits
        )
      :
        pPFParticle(input_PFParticle),
        pMetaData(input_MetaData),
        pVertex(input_Vertex),
        pClusters(input_Clusters)
    {
      pPFParticleID = pPFParticle->Self();
      pPdgCode = pPFParticle->PdgCode();

      //Get recob::Shower/Track
      const unsigned int nShowers(input_Showers.size());
      const unsigned int nTracks(input_Tracks.size());
      if(nShowers+nTracks != 1) mf::LogDebug("SinglePhoton") << "  No just one shower/track is associated to PFParticle " << pPFParticleID << "\n";

      pHasShower = nShowers;
      pHasTrack = nTracks;
      if(pHasShower == 1) pShower=input_Showers.front();
      if(pHasTrack == 1) pTrack=input_Tracks.front();


      //fill the vertex info.
      if (!pVertex.empty()){
        const art::Ptr<recob::Vertex> vertex = * (pVertex.begin());
        vertex->XYZ(pVertex_pos);
      } 

      //fill pPFPHits from each clusters;
      for(size_t index=0; index< pClusters.size(); ++index){
        auto cluster = pClusters[index];
        std::vector< art::Ptr<recob::Hit  > > temp_hits = input_Hits.at(cluster.key());
        pPFPHits.insert(pPFPHits.end(), temp_hits.begin(), temp_hits.end()); 
      }

      //get ancestor for a 1st generation PFParticle
      if(pPFParticle->IsPrimary()){ 
        pAncestor = pPFParticle;
        pAncestorID = -1;
      }

      //fill in some booleans
      for(auto &meta: pMetaData){
        std::map<std::string, float> propertiesmap  = meta->GetPropertiesMap();
        if(propertiesmap.count("NuScore")==1) pNuScore = propertiesmap["NuScore"];
        if(propertiesmap.count("TrackScore")==1) pTrackScore = propertiesmap["TrackScore"];
        if(propertiesmap.count("IsNeutrino")==1) pIsNeutrino = true; 
        if(propertiesmap.count("IsClearCosmic")==1) pIsClearCosmic = true; 
      }

    };

    art::Ptr< recob::PFParticle > pPFParticle;//d
    art::Ptr< recob::Shower>  pShower;//d with 0 or 1 element
    art::Ptr< recob::Track >  pTrack;//d with 0 or 1 element
    art::Ptr< recob::Slice >  pSlice;//d in helper_connector.h get the id from pSlice->key()
    art::Ptr< recob::PFParticle > pAncestor;//d found by tracing Parent()

    art::Ptr<anab::ParticleID>  pParticleID;//d for track only;
    art::Ptr< simb::MCTruth >  pMCTruth;

    std::vector< art::Ptr< larpandoraobj::PFParticleMetadata > > pMetaData;//d
    std::vector< art::Ptr< recob::Vertex > > pVertex;//d
    std::vector< art::Ptr< recob::Hit > >  pSliceHits;//d
    std::vector< art::Ptr< recob::Hit > >  pPFPHits;//d
    std::vector< art::Ptr< recob::Cluster > > pClusters;//d
    std::vector<art::Ptr<anab::Calorimetry>> pCalorimetries;//d
    std::vector< art::Ptr< recob::SpacePoint > > pSpacePoints;
    std::vector< art::Ptr< simb::MCParticle > > pMCParticles;


    //set methods
    void set_NuScore (const double input_score){ pNuScore = input_score; }
    //      void set_TrackScore (const double input_score){ pTrackScore = input_score; }

    //      void set_IsNeutrino (const bool input_bool){ pIsNeutrino = input_bool; }
    //      void set_IsClearCosmic (const bool input_bool){ pIsClearCosmic = input_bool; }
    void set_IsNuSlice (const bool input_bool){ pIsNuSlice = input_bool; }
    void set_HasPID (const bool input_bool){ pHasPID = input_bool; }

    //      void set_HasShower (const int input_number){ pHasShower = input_number; }
    //      void set_HasTrack (const int input_number){ pHasTrack = input_number; }
    //      void set_PdgCode (const int input_number){ pPdgCode = input_number; }
    //      void set_PFParticleID (const int input_number){ pPFParticleID = input_number; }
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
  void PPFP_FindSliceIDandHits(std::vector< PandoraPFParticle > & PPFPs, art::Ptr< recob::Slice >  slice, const std::vector<art::Ptr<recob::PFParticle> > PFP_in_slice, const std::vector<art::Ptr<recob::Hit> > Hit_inslice){

    int pfp_size = PPFPs.size();
    for( auto pfp : PFP_in_slice){
      for(int index = 0; index < pfp_size; index++){
        if(PPFPs[index].get_SliceID() > -1 ) continue;//slice# is found already
        if( (PPFPs[index].pPFParticle)->Self() == pfp->Self() ){
          PPFPs[index].pSlice = slice;
          PPFPs[index].pSliceHits = Hit_inslice;
          PPFPs[index].set_SliceID( slice.key() );
          break;
        }
      }
    }
  }

  //refill pNuScore and pIsNuSlice
  int DefineNuSlice(std::vector< PandoraPFParticle > & PPFPs){

    int pfp_size = PPFPs.size();
    double best_nuscore = 0;
    int best_nuscore_SliceID = 0;
    std::vector< int > IDs;

    for(int index = 0; index < pfp_size; index++){
      PandoraPFParticle* temp_p = &PPFPs[index];
      if(temp_p->get_IsNeutrino()){
        int temp_sliceID = temp_p->get_SliceID();
        //add one if not found;
        if(!std::count(IDs.begin(), IDs.end(), temp_sliceID) ) IDs.push_back(temp_sliceID);
        if(best_nuscore < temp_p->get_NuScore() ){
          best_nuscore = temp_p->get_NuScore();

          best_nuscore_SliceID = temp_p->get_SliceID();

        }
      }
    }

    //now markdown all pfparticles in slice
    //re-set pNuScore and pIsNuSlice
    for(int index = 0; index < pfp_size; index++){
      PandoraPFParticle* ppfp = &PPFPs[index];
      if( std::count(IDs.begin(), IDs.end(), ppfp->get_SliceID()) ) ppfp->set_IsNuSlice(true);
    }

    return best_nuscore_SliceID;
  }



//2. Trackers to find the correct PandoraPFParticle
  PandoraPFParticle *PPFP_GetPPFPFromShower( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Shower> pShower){
    int pfp_size = PPFPs.size();
    for(int index = 0; index < pfp_size; index++){
      if(PPFPs[index].get_HasShower() != 1 ) continue;
//      std::cout<<"CHECK Shower start "<<pShower->ShowerStart().X()<<" vs "<<PPFPs[index].pShower->ShowerStart().X()<<std::endl;
      //CHECK, this works, but there maybe a better way;
      if((pShower->ShowerStart() == PPFPs[index].pShower->ShowerStart())
      && (pShower->Direction() == PPFPs[index].pShower->Direction())){
        return &PPFPs[index];
      }
    }
    std::cout<<"Error, no PFParticle matched to shower, returning the first element"<<std::endl;
    return &PPFPs[0];
  }

  PandoraPFParticle *PPFP_GetPPFPFromTrack( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Track> pTrack){
    int pfp_size = PPFPs.size();
    for(int index = 0; index < pfp_size; index++){
      if(PPFPs[index].get_HasTrack() != 1 ) continue;
      if((pTrack->StartDirection() == PPFPs[index].pTrack->StartDirection())
      && (pTrack->EndDirection() == PPFPs[index].pTrack->EndDirection())){
        return &PPFPs[index];
      }
    }
    std::cout<<"Error, no PFParticle matched to track, returning the first element"<<std::endl;
    return &PPFPs[0];
  }

  PandoraPFParticle *PPFP_GetPPFPFromPFID( std::vector< PandoraPFParticle > & PPFPs, int id){
    int pfp_size = PPFPs.size();
    for(int index = 0; index < pfp_size; index++){
      if(PPFPs[index].get_PFParticleID() == id ){
        return &PPFPs[index];
      }
    }
    std::cout<<"Error, no PFParticle matched to track, returning the first element"<<std::endl;
    return &PPFPs[0];
  }





}
