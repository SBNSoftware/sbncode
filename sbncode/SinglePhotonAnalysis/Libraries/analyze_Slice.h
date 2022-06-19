namespace single_photon
{
    //called once per event to get all the slice info
    //fills a map between the neutrino score for the slice and the primary reco PFP
    //loops over all PFP's to find the primary and then associate to a slice
    void SinglePhoton::AnalyzeSlices(std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap,
            PFParticleIdMap &pfParticleMap,
            std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> &primaryPFPSliceIdVec,
            std::map<int, double> &sliceIdToNuScoreMap,
            std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
            std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap,
            std::map<art::Ptr<recob::PFParticle>,bool>& PFPToNuSliceMap,
            std::map<art::Ptr<recob::PFParticle>,double>& PFPToTrackScoreMap){


        //std::vector<std::pair<art::Ptr<recob::PFParticle>, int>> primaryPFPSliceIdVec; //maps a primary PFP to a slice index
        // std::map<int, double> sliceIdToNuScoreMap; //maps a slice index to the associated neutrino score
        std::vector<art::Ptr<recob::PFParticle>> clearCosmicPFP;
        std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind


        std::vector<double> nuscore_slices; //this is a temporary vector to store neutrino score per slice for this event
        std::map<int, bool> sliceIdToNuSliceMap; //this is a temporary vector to store neutrino score per slice for this event
        //std::vector<art::Ptr<recob::PFParticle>> primary_pfps; //store the primary PFP for each slice
        // sliceIdToPFPMap.clear(); //clear between events

        /*
         * Grabbed this info from Giuseppe:
         * Here's the structure of these Metadata
         * Primary PfParticles are either
         * 1) IsClearCosmic = 1 (for unambiguous cosmics)
         * 2) NuScore = 0.108586, SliceIndex = 1 (for cosmic slices)
         * 3) IsNeutrino = 1, NuScore = 0.170914, SliceIndex = 2 (for the nu slice)
         * Then, for PfParticles that are daughter of the nu, the track score is saved, e.g.:
         * 4) TrackScore = 0.671488
         * PfParticles that are not primary and that are not daugthers of the neutrino have empty Metadata
         */


        //for all PFP metadata in the event
        for (auto pair: pfParticleToMetadataMap){
            std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatalist= pair.second; //get the metadata
            art::Ptr<recob::PFParticle> pfp = pair.first; //get the corresponding PFP

            //will be empty in circumstances outlined above, not every PFP has stored metadata
            if (!metadatalist.empty()){

                //std::cout<<"metadatalist not empty for pfp with index and pdg code: "<<pfp->Self()<<"/"<<pfp->PdgCode()<<", primary = "<<pfp->IsPrimary()<<std::endl;

                //for each PFP per event
                for(art::Ptr<larpandoraobj::PFParticleMetadata> data:metadatalist){

                    //const pandora::PropertiesMap &pfParticlePropertiesMap(metadata->GetPropertiesMap()); 

                    //get the metadata properties
                    std::map<std::string, float> propertiesmap  = data->GetPropertiesMap();
                    //std::cout<<"the number of items in the metadata properties map is "<<propertiesmap.size()<<std::endl;
                    int temp_ind = -1;
                    double temp_score = -1.0;
                    int clear_cosmic = -1;
                    bool is_nuslice = false;
                    //for each of the things in the list
                    for (auto it:propertiesmap ){
                        // std::cout << "  - " << it.first << " = " << it.second << std::endl;
                        if (it.first == "SliceIndex"){
                            temp_ind = it.second;
                            // std::cout << "  - " << it.first << " = " << it.second << std::endl;
                        }
                        //store the neutrino score for each slice
                        if (it.first == "NuScore"){
                            nuscore_slices.push_back(it.second);
                            temp_score = it.second;
                            //std::cout << "  - " << it.first << " = " << it.second << std::endl;
                            //if it's the neutrino score it also means it's the primary particle so save the pfp
                            //primary_pfps.push_back(pfp);
                            //pfParticleToNuScoreMap[pfp] = it.second;


                        }
                        if (it.first == "IsClearCosmic"){
                            clear_cosmic = 1;
                        }
                        if(it.first == "IsNeutrino"){
                            is_nuslice = true;
                        }
                        if(it.first == "TrackScore"){
                            //std::cout << "  - " << it.first << " = " << it.second << std::endl;
                            //std::cout << "  - "<<pfp->Self()<<" "<< it.first << " = " << it.second << std::endl;
                            PFPToTrackScoreMap[pfp] = it.second;
                        }

                    }//for each item in properties map

                    //if there is a neutrino score it's the primary PFP, so save the score+slice info
                    if(temp_score != -1.0){
                        primaryPFPSliceIdVec.push_back(std::pair(pfp, temp_ind));
                        //primaryToSliceIdMap[pfp] = temp_ind;
                        sliceIdToNuScoreMap[temp_ind] = temp_score;
                        if(m_is_verbose)std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t found primary PFP at index "<<pfp->Self()<<" at slice Index "<<temp_ind<<" has a nu score of"<<temp_score<<std::endl;

                    }
                    if( clear_cosmic> 0){
                        //std::cout<<"clear cosmic with PFP id "<<pfp->Self()<<std::endl;
                        clearCosmicPFP.push_back(pfp);
                        PFPToClearCosmicMap[pfp] = true;
                    }else{
                        PFPToClearCosmicMap[pfp] = false;

                    }

                    sliceIdToNuSliceMap[temp_ind] = is_nuslice;

                }//for each PFP/metadata

            }//if the list isn't empty

        }//for all PFP/metadata in event

        //now we have all the primary pfp's and the corresponding slices+scores
        //the next step is to look at all the pfp's in the event, find the primary, and then store the slice ind
        // std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind

        //for all pfp's in the event
        //std::cout<<"looking at all PFP's"<<std::endl;
        //std::map<int, std::vector<art::Ptr<recob::PFParticle>>> sliceIdToPFPMap;

        //for (unsigned int i = 0; i< pfParticleVector.size(); i++){}
        for(auto item: pfParticleMap){
            art::Ptr<recob::PFParticle> start_pfp = item.second;
            //no longer skipping pfp's that are clear cosmics
            //auto iter = find(clearCosmicPFP.begin(), clearCosmicPFP.end(), start_pfp);
            //if(iter != clearCosmicPFP.end()) continue;

            // std::cout<<"START: looking for match for pfp - track id/pdg code "<<start_pfp->Self()<<"/"<<start_pfp->PdgCode()<<std::endl; 

            art::Ptr<recob::PFParticle> this_pfp = start_pfp;
            art::Ptr<recob::PFParticle> parent_pfp ;

            //if this is a primary particle, skip next part
            if(!this_pfp->IsPrimary()){
                parent_pfp = pfParticleMap[this_pfp->Parent()];
                //std::cout<<"parent of start particle is "<<parent_pfp->Self()<<"/"<<parent_pfp->PdgCode()<<std::endl;   
                //if not primary, iterate up parent chain
                while(!parent_pfp->IsPrimary()){
                    //  std::cout<<"not primary - track id/pdg code "<<parent_pfp->Self()<<"/"<<parent_pfp->PdgCode()<<std::endl; 
                    //std::cout<<"iterating, current particle has index "<<parent_pfp->Self()<<std::endl;
                    parent_pfp = pfParticleMap[this_pfp->Parent()];
                    this_pfp = parent_pfp;
                }//while not primary, iterate up chain
            } else{
                parent_pfp = start_pfp;
            }

            //std::cout<<"this particle was primary at track id/pdg code "<<parent_pfp->Self()<<"/"<<parent_pfp->PdgCode()<<std::endl; 

            //get slice id for this primary
            int slice_id = -1;//initialize to invalid

            //for each primary pfp
            for(unsigned int j = 0; j <primaryPFPSliceIdVec.size(); j++){
                //if the parent primary matches to a primary in the list, save the slice id
                if (parent_pfp == primaryPFPSliceIdVec[j].first){
                    slice_id = primaryPFPSliceIdVec[j].second;
                    break;
                }
            }//for all primary PFP's in event

            //store original pfp and it's slice id
            if(slice_id < 0 ){
//                if(m_is_verbose)std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t no matching slice found for this PFP with primary (parent) id "<<parent_pfp->Self()<<std::endl;
            } else {
                // allPFPSliceIdVec[i] = std::pair(start_pfp,slice_id);
                allPFPSliceIdVec.push_back(std::pair(start_pfp,slice_id));
                if(m_is_verbose)std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t PFParticle "<<start_pfp->Self()<<"pdg: ("<<start_pfp->PdgCode()<<") at slice "<<slice_id<<" has PFP with primary id "<<parent_pfp->Self()<<std::endl;
                // PFPToSliceIdMap[start_pfp] = slice_id; 
            }
			//std::cout<<"CHECK "<<__LINE__<<" work on slice ID: "<<slice_id<<std::endl;

            PFPToSliceIdMap[start_pfp] = slice_id;
            //std::cout<<"storing PFP with id "<<start_pfp->Self()<<" and slice id "<<slice_id<<" in  PFPToSliceIdMap"<<std::endl; 
            PFPToNuSliceMap[start_pfp] = sliceIdToNuSliceMap[slice_id];


            // sliceIdToPFPMap[slice_id].push_back(start_pfp);
        }//for all pfp's in the event

//CHECK_not used        for (auto pair: sliceIdToPFPMap){
//             std::cout<<"in slice ID "<<pair.first<<" there are "<<pair.second.size()<<" PFP's"<<std::endl;
//        } 

//CHECK , useless?
//for (auto pair:PFPToNuSliceMap ){
//            if(pair.second == true){
//                std::cout<<"stored in PFPToNuSliceMap for pfp "<<pair.first->Self()<<", isNeutrino = "<<pair.second<<std::endl;
//            }
//        }

        /*
         * store stuff in the output tree
         */

        m_reco_slice_num = nuscore_slices.size()+1;//the number of slices also corresponds to the number of neutrino scores

        //currently this is junk, just a placeholder
        //std::cout<<"saving the info for "<<m_reco_slice_num<<" slices"<<std::endl;
        this->ResizeSlices(m_reco_slice_num); 
        m_reco_slice_nuscore = nuscore_slices;

        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of clear cosmic PFP's in the event is "<<clearCosmicPFP.size()<<std::endl;
        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of items in the primary to slice vector is "<<primaryPFPSliceIdVec.size()<<" and in the slice to score map is "<<sliceIdToNuScoreMap.size()<<std::endl;
        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t the number of PFP's matched to a slice is "<< allPFPSliceIdVec.size()<<"/"<< pfParticleMap.size()<<std::endl;
        if ((clearCosmicPFP.size() +  allPFPSliceIdVec.size())!= pfParticleMap.size()){
            std::cout<<"BIG ERROR, UNACCOUNTED FOR PFP's, (clearCosmicPFP.size() +  allPFPSliceIdVec.size())!= pfParticleMap.size())"<<std::endl;
        }

        if (PFPToSliceIdMap.size() != pfParticleMap.size()){
            std::cout<<"BIG ERROR, PFPToSliceIdMap.size() != pfParticleMap.size(): "<<std::endl;
            std::cout<<PFPToSliceIdMap.size()<<"/"<<pfParticleMap.size()<<std::endl;
        }

        if (PFPToClearCosmicMap.size() != pfParticleToMetadataMap.size()){
            std::cout<<"BIG ERROR, PFPToClearCosmicMap.size() != pfParticleToMetadataMap.size(): "<<std::endl;
            std::cout<<PFPToClearCosmicMap.size()<<"/"<<pfParticleToMetadataMap.size()<<std::endl;

        } 

        //std::vector<std::vector<int>> pfp_pdg_slice; 
        //for(auto item: allPFPSliceIdVec){
        //std::cout<<"the pfp with id "<<item.first->Self()<<" is associated to slice "<<item.second<<std::endl;
        //      pfp_pdg_slice[item.second].push_back(item.first->PdgCode());
        // }
    }


    std::vector<int>  SinglePhoton::GetPFPsPerSlice( std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap){
        std::vector<int> sliceIdToNumPFPsvec(m_reco_slice_num, 0);

        //std::cout<<"starting to look at the PFP's per slice"<<std::endl;

        //if the map isn't filled, return 0 PFP's per slice
        if(PFPToSliceIdMap.size() < 1){ 
            std::cout<<"Error, no PFP's in PFPToSliceIdMap, size = "<<PFPToSliceIdMap.size()<<std::endl;
            return sliceIdToNumPFPsvec;
        }

        //for all PFP's
		//CHECK the following block crashes the code
		//FIXED by fixing m_reco_slice_num (add 1)
        for (auto pair:PFPToSliceIdMap ){

            //get the slice and increment the vector
            int slice_id = pair.second;
			if(slice_id > -1) sliceIdToNumPFPsvec[slice_id]++;
//			std::cout<<"CHECK found 1 pfp particle at slice id "<<slice_id<<"/"<<m_reco_slice_num<<std::endl;
        }

        return sliceIdToNumPFPsvec;

    }

    std::vector<int> SinglePhoton::GetNumShowersPerSlice(std::map< art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& showerToPFParticleMap, 
            std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap){
        std::vector<int> sliceIdToNumShowersvec(m_reco_slice_num, 0);

        /*
           int max_pfp = -1;
           for(auto it = PFPToSliceIdMap.begin(); it != PFPToSliceIdMap.end(); ++it ) {
           int this_pfp = it->first->Self();
           if (this_pfp> max_pfp) {
           max_pfp = this_pfp;
           }

           }
           std::cout<<"the max PFP id in PFPToSliceIdMap is  = "<< max_pfp<<std::endl;
           */
        int cosmics = 0;
        //for each shower
        for (auto pair: showerToPFParticleMap){
            art::Ptr<recob::PFParticle> pfp = pair.second;

            //check if this is a valid PFP
            //if (pfp->self() > PFPToSliceIdMap.size()){
            // std::cout<<"ERROR, this pfp is out of bounds "

            // }

            //find slice corresponding to PFP
            //have to check if it's in the map otherwise it returns 0 which is misleading
            if (PFPToSliceIdMap.find(pfp) != PFPToSliceIdMap.end()){
                //get the slice id

                int slice_id = PFPToSliceIdMap[pfp];
                if (slice_id > -1){
                    //std::cout<<"looking at shower id  "<<pair.first->ID()<<" with pfp "<<pair.second->Self()<<" in slice "<<slice_id<<std::endl;
                    //incrmement number of tracks per slice
                    int num =  sliceIdToNumShowersvec[slice_id];
                    sliceIdToNumShowersvec[slice_id]= ++num;
                }//if neutrino slice id
                else{
                    cosmics++;
                }
            }//if there's a fpf-slice match
        }//for each shower

        unsigned int count = 0;
        for (auto numshowers: sliceIdToNumShowersvec){
            count += numshowers;

        }
        return sliceIdToNumShowersvec;
    }

	std::vector<int> SinglePhoton::GetNumTracksPerSlice(
			std::map< art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>>& trackToPFParticleMap, 
			std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap){
        std::vector<int> sliceIdToNumTracksvec(m_reco_slice_num, 0);

        // std::cout<<"looking at number of tracks per slice"<<std::endl;
        int cosmics = 0;
        //for each track
        for (auto pair: trackToPFParticleMap){
            art::Ptr<recob::PFParticle> pfp = pair.second;
            //find slice corresponding to PFP
            //have to check if it's in the map otherwise it returns 0 which is misleading
            if (PFPToSliceIdMap.find(pfp) != PFPToSliceIdMap.end()){
                //get the slice id

                int slice_id = PFPToSliceIdMap[pfp];
                if (slice_id > -1){
                    // std::cout<<"looking at track id  "<<pair.first->ID()<<" with pfp "<<pair.second->Self()<<" in slice "<<slice_id<<std::endl;
                    //incrmement number of tracks per slice
                    int num =  sliceIdToNumTracksvec[slice_id];
                    sliceIdToNumTracksvec[slice_id]= ++num;
                }//if neutrino slice id
                else{
                    cosmics++;
                }
            }//if there's a fpf-slice match
        }//for each track
        unsigned int count = 0;
        for (auto numtracks: sliceIdToNumTracksvec){
            count += numtracks;

        }
        if (count != (trackToPFParticleMap.size()-cosmics)){
            std::cout<<"Error, number of showers in sliceIdToNumTracksvec is "<<count<<" but trackToPFParticleMap.size() for neutrino showers = "<<trackToPFParticleMap.size()-cosmics<<std::endl;
        }



        return sliceIdToNumTracksvec;
    }



    //if 1g1p
    //is there at least 1 reco track and shower
    //if there yes, are they in the same slice
    //if yes, what is the slice score?
    //if yes, was it a cosmic slice?
    //if yes is that the neutrino slice?
    //if they are in different slices
    //was the shower/track:
    //in the neutrino slice?
    //in a cosmic slice?

    //if 1g0p 
    //if there is 1 shower
    //what is the slice id? slice score? cosmic?
    //if missing, not-recoed         

    //if missing atleast 1shower, not recoed


    //for a given signal def (`ncdelta` for now) , finds the MCParticles in event
    //loops over association between reco tracks/showers to get associated slice(s)
    //can also look at things like shower energy, conversion length, etc.
    void SinglePhoton::AnalyzeRecoMCSlices(std::string signal_def, 
			std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,
            std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle> > & showerToPFParticleMap, 
            std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec, 
            std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,
            std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle> > & trackToNuPFParticleMap,
            std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap,
            std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap){

        m_reco_slice_num_pfps = GetPFPsPerSlice(PFPToSliceIdMap); //the total number of PFP's per slice
		m_reco_slice_num_showers = GetNumShowersPerSlice(showerToPFParticleMap, PFPToSliceIdMap); //the subset of PFP's that are showers
        m_reco_slice_num_tracks = GetNumTracksPerSlice(trackToNuPFParticleMap, PFPToSliceIdMap);


        //first check if in the event there's a match to a given signal
        if(signal_def == "ncdelta"){
			//@para updated in the AnalyzeMCTruths function @ analyze_MCTruth.h
                std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t looking for signal def "<<signal_def<<", m_mctruth_is_delta_radiative = "<<m_mctruth_is_delta_radiative<<std::endl; 

                std::vector<int> matched_shower_ids;
                //first look for sim showers
                for (unsigned int j = 0; j< m_sim_shower_parent_pdg.size(); j++){
                    int parent= m_sim_shower_parent_pdg[j];
                    int pdg =  m_sim_shower_pdg[j];

                    //if this sim shower is a photon and it's primary (parent pdg is -1)
                    if(parent == -1 && pdg ==22){
                        //first check that this particle isn't alread saved
                        //use map from track ID to get MCP
                        //if this shower is matched to a recob:shower
                        if (m_sim_shower_matched[j] > 0 &&  m_reco_shower_energy_max[j] >20){
                            int matched_shower_id = m_sim_shower_trackID[j];


                            //if this shower isn't already stored
                            if (std::find(matched_shower_ids.begin(), matched_shower_ids.end(), matched_shower_id) == matched_shower_ids.end()){
                                matched_shower_ids.push_back(matched_shower_id);
                                m_matched_signal_shower_overlay_fraction.push_back(m_sim_shower_overlay_fraction[j]);
                                //m_matched_signal_shower_conversion_length;
                                m_matched_signal_shower_true_E.push_back(m_sim_shower_energy[j]);
                                m_matched_signal_shower_nuscore.push_back( m_reco_shower_nuscore[j]);
                                int id = m_reco_shower_sliceId[j];
                                std::cout<<"found matched photon shower in slice "<<id<<" with m_sim_shower_energy[j] = "<<m_sim_shower_energy[j]<<std::endl;

                                m_matched_signal_shower_sliceId.push_back(id);


                                m_matched_signal_shower_is_clearcosmic.push_back( m_reco_shower_isclearcosmic[j]);
                                m_matched_signal_shower_is_nuslice.push_back(m_reco_shower_is_nuslice[j]);
				// num track/shower in slice here is in reverse order
                                m_matched_signal_shower_tracks_in_slice.push_back(m_reco_slice_num_showers[id]);
                                m_matched_signal_shower_showers_in_slice.push_back(m_reco_slice_num_tracks[id]);
                                // std::cout<<"found signal photon shower pdg"<< m_sim_shower_pdg[j]<<"and is in neutrino slice =  "<< m_sim_shower_is_nuslice[j]<<std::endl;

                            }//if not already stored
                        }//if matched to a reco shower >20MeV
                    }//if it's a photon from the neutrino interaction
                }//for all sim showers

                m_matched_signal_shower_num = m_matched_signal_shower_true_E.size();

				//NEXT, same procedure for tracks
                std::vector<int> matched_track_ids;
                for (unsigned int k = 0; k< m_sim_track_parent_pdg.size(); k++){
                    int parent= m_sim_track_parent_pdg[k];
                    int pdg =  m_sim_track_pdg[k];

                    int matched_track_id = m_sim_track_trackID[k];

                    //if this sim track is a proton and it's primary (parent pdg is -1)
                    if((parent == -1 ||parent == 12 || parent ==14 ) && pdg == 2212){

                        if (m_sim_track_matched[k] > 0){

                            if (std::find(matched_track_ids.begin(), matched_track_ids.end(), matched_track_id) == matched_track_ids.end()){
                                matched_track_ids.push_back(matched_track_id);


                                // m_matched_signal_track_overlay_fraction.push_back(m_sim_track_overlay_fraction[j]);
                                m_matched_signal_track_true_E.push_back(m_sim_track_energy[k]);
                                m_matched_signal_track_nuscore.push_back( m_reco_track_nuscore[k]);
                                m_matched_signal_track_sliceId.push_back(m_reco_track_sliceId[k]);
                                m_matched_signal_track_is_clearcosmic.push_back( m_reco_track_isclearcosmic[k]);
                                m_matched_signal_track_is_nuslice.push_back(m_reco_track_is_nuslice[k]);
                               
                                int id = m_reco_track_sliceId[k];
                                m_matched_signal_track_tracks_in_slice.push_back(m_reco_slice_num_tracks[ id]);
                                m_matched_signal_track_showers_in_slice.push_back(m_reco_slice_num_showers[ id]);

                            }

                        }//if matched
                    }//if proton from neutrino interaction
                }//for all sim tracks
                m_matched_signal_track_num = m_matched_signal_track_true_E.size();
        }//end of "ncdelta" scenario
  	
  	//brief summary
        if (m_matched_signal_shower_num > 1) m_multiple_matched_showers = true;
        if (m_matched_signal_track_num > 1) m_multiple_matched_tracks = true;
        if (m_matched_signal_shower_num == 0)  m_no_matched_showers = true;

        //check if either 1g1p or 1g0p topology
        if (m_matched_signal_shower_num ==1  && m_matched_signal_track_num ==1){//1g1p
            //check if same slice
                m_is_matched_1g1p = true;
            if ( m_matched_signal_track_sliceId[0] == m_matched_signal_shower_sliceId[0]){
                m_reco_1g1p_is_same_slice = true;
                m_reco_1g1p_is_nuslice = m_matched_signal_shower_is_nuslice[0];
                m_reco_1g1p_nuscore = m_matched_signal_shower_nuscore[0];
            } else{
                m_reco_1g1p_is_multiple_slices = true;
            }
        }else if(m_matched_signal_shower_num ==1  && m_matched_signal_track_num ==0){//1g0p
            m_reco_1g0p_is_nuslice = m_matched_signal_shower_is_nuslice[0];
            m_reco_1g0p_nuscore =  m_matched_signal_shower_nuscore[0];
            m_is_matched_1g0p = true;

        }

  	//CHECK, text reports
		if(m_is_verbose){
			if (m_reco_1g1p_is_nuslice){
				int nu_id = m_matched_signal_shower_sliceId[0];
				for (auto pair: PFPToSliceIdMap){
					auto pfp = pair.first;
					auto id = pair.second;
					if (id == nu_id){
						std::cout<<"the pfp in this nu slice with id "<<id <<" is "<<pfp->Self()<<" with pdg "<<pfp->PdgCode()<<std::endl;
					}
				}
			}

			if (m_reco_1g1p_is_nuslice && ( m_matched_signal_shower_tracks_in_slice[0]>1 ||  m_matched_signal_shower_showers_in_slice[0]>1) ){
				std::cout<<"found reco 1g1p with "<<  m_matched_signal_shower_showers_in_slice[0]<<" showers and "<< m_matched_signal_shower_tracks_in_slice[0]<<"tracks in the same slice in run = "<<m_run_number <<", subrun = "<<m_subrun_number  <<", event = "<<m_event_number<<std::endl;

			}

			for(unsigned int i= 0; i<m_reco_slice_num_pfps.size(); i++){
				std::cout<<"The number of PFP's in slice: " << i<<std::endl;

				std::cout<< "-- is "<<m_reco_slice_num_pfps[i]<<std::endl;
				std::cout<<"-- -- of which "<< m_reco_slice_num_showers[i]<< " are showers and "<< m_reco_slice_num_tracks[i] <<" are tracks"<<std::endl;
				if (m_reco_slice_num_showers[i] + m_reco_slice_num_tracks[i] != m_reco_slice_num_pfps[i]){
					std::cout<<"ERROR, mismatching numbers of PFPs"<<std::endl;
				}
			}
		}
        std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t "<<__FUNCTION__<<" is done "<<std::endl;

    }//findslice

}
