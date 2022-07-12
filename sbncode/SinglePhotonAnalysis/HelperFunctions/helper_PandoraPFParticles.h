namespace single_photon
{//helper functions exclusively for PandoraPFParticles
	
//1. Filler for each PandoraPFParticle
	//find pAncestor & pAncesotrID
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
				//			std::cout<<PPFPs[jndex].pPFParticleID<<" Trace up a generation parent "<<temp_parent_id<<std::endl;

			}
		}
	}

	//Fill pSlice, pHits, pSliceID
	void PPFP_FindSliceIDandHits(std::vector< PandoraPFParticle > & PPFPs, art::Ptr< recob::Slice >  slice, const std::vector<art::Ptr<recob::PFParticle> > PFP_in_slice, const std::vector<art::Ptr<recob::Hit> > Hit_inslice){

		int pfp_size = PPFPs.size();
		for( auto pfp : PFP_in_slice){
			for(int index = 0; index < pfp_size; index++){
				//std::cout<<"CHECK slice match"<<(PPFPs[index].pPFParticle)->Self()<< " and "<<pfp->Self()<<std::endl;
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
//		std::map< int, double > temp_ID_nuscore_map;
		double best_nuscore = 0;
		int best_nuscore_SliceID = 0;
		std::vector< int > IDs;

		for(int index = 0; index < pfp_size; index++){
			PandoraPFParticle* temp_p = &PPFPs[index];
			if(temp_p->get_IsNeutrino()){
				int temp_sliceID = temp_p->get_SliceID();
				//add one if not found;
				if(!std::count(IDs.begin(), IDs.end(), temp_sliceID) ) IDs.push_back(temp_sliceID);
//				temp_ID_nuscore_map[ temp_p.get_SliceID()] = temp_p.get_NuScore();
				if(best_nuscore < temp_p->get_NuScore() ){
					best_nuscore = temp_p->get_NuScore();

					best_nuscore_SliceID = temp_p->get_SliceID();

////					std::cout<<__FUNCTION__<<"CHECK Best nu score "<< best_nuscore<<" at slice "<<temp_p.get_SliceID()<<std::endl;
				}
			}
		}

		//now markdown all pfparticles in slice
		//re-set pNuScore and pIsNuSlice
		for(int index = 0; index < pfp_size; index++){
			PandoraPFParticle* ppfp = &PPFPs[index];
			if( std::count(IDs.begin(), IDs.end(), ppfp->get_SliceID()) ) ppfp->set_IsNuSlice(true);
//			if( PPFPs[index].get_SliceID() == best_nuscore_SliceID){
//				PPFPs[index].set_IsNuSlice( true );
//				PPFPs[index].set_NuScore( best_nuscore );//CHECK over-write the original score, if there is any;
//				std::cout<<__FUNCTION__<<" CHECK Set nu slice "<<PPFPs[index].get_SliceID()<<" with score "<<PPFPs[index].get_NuScore()<<std::endl;
//			}//CHECK, same event, different slices, still different score.
		}

		return best_nuscore_SliceID;
	}



//2. Trackers to find the correct PandoraPFParticle
	PandoraPFParticle *PPFP_GetPPFPFromShower( std::vector< PandoraPFParticle > & PPFPs, art::Ptr<recob::Shower> pShower){
		int pfp_size = PPFPs.size();
		for(int index = 0; index < pfp_size; index++){
			if(PPFPs[index].get_HasShower() != 1 ) continue;
//			std::cout<<"CHECK Shower start "<<pShower->ShowerStart().X()<<" vs "<<PPFPs[index].pShower->ShowerStart().X()<<std::endl;
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



}
