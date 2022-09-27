#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"

namespace single_photon
{
  PandoraPFParticle::PandoraPFParticle( 
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

  }

  //helper functions exclusively for PandoraPFParticles
  void PPFP_FindAncestor ( std::vector< PandoraPFParticle > & PPFPs){

    std::map< size_t, art::Ptr<recob::PFParticle>> pfParticleMap;
    int pfp_size = PPFPs.size();
    //build ID-PFParticle map;
    for(int index = 0; index < pfp_size; index++){
      PandoraPFParticle temp_pfp = PPFPs[index];
      if (!pfParticleMap.insert(std::map< size_t, art::Ptr<recob::PFParticle>>::value_type((temp_pfp.pPFParticle)->Self(), temp_pfp.pPFParticle)).second){
        throw cet::exception("SinglePhoton") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
      }
    }

    //trace up parents
    for(int jndex = 0; jndex < pfp_size; jndex++){
      art::Ptr< recob::PFParticle > temp_pfp = PPFPs[jndex].pPFParticle;

      PPFPs[jndex].set_AncestorID(temp_pfp->Self() );
      PPFPs[jndex].pAncestor  = pfParticleMap[temp_pfp->Self()];
      if(temp_pfp->IsPrimary()) continue;//skip PFP without parent is a parent of itself

      while(!(PPFPs[jndex].pAncestor)->IsPrimary()){//1+ gen. parent

        int temp_parent_id = PPFPs[jndex].pAncestor->Parent();
        PPFPs[jndex].set_AncestorID( temp_parent_id );
        PPFPs[jndex].pAncestor  = pfParticleMap[temp_parent_id];
        //      std::cout<<PPFPs[jndex].pPFParticleID<<" Trace up a generation parent "<<temp_parent_id<<std::endl;

      }
    }

  }


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
  PandoraPFParticle* PPFP_GetPPFPFromShower( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Shower> pShower){
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


  PandoraPFParticle* PPFP_GetPPFPFromTrack( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Track> pTrack){
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


  PandoraPFParticle* PPFP_GetPPFPFromPFID( std::vector< PandoraPFParticle > & PPFPs, int id){
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
