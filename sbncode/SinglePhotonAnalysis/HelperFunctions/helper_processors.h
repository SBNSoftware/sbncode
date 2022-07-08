namespace single_photon
{
	//--------------------- a template to get vectors from labels -------------------------
	//1. Scratch from Handle, and return a equivalently useful vector.
	//  Sample usage:
	//		art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); //yea, this should be gone;
	//		recob::Hit dummy_hit;//This is to specify the template;
	//		std::vector<art::Ptr<recob::Hit>> hitVector = VectorFromLabel(dummy_hit, evt, m_hitfinderLabel);

//	template <typename recob_object>//A helper template that allows you to make compliated types.
//		struct temporary_types{ 
//			using type1 = std::vector<art::Ptr<recob_object>>;
//			using type2 = art::ValidHandle<std::vector<recob_object>>;
//			using type3 = std::vector<recob_object>;
//		};
//
//	template <class recob_object>//ref_type is only used to identify the temporary_types late.
//		typename temporary_types<recob_object>::type1 VectorFromLabel(recob_object ref_type, const art::Event &evt, std::string &label){
//
//			typename temporary_types<recob_object>::type2 const & Handle = evt.getValidHandle<typename temporary_types<recob_object>::type3>(label);
//			typename temporary_types<recob_object>::type1 Vector;
//			art::fill_ptr_vector(Vector,Handle);
//			return Vector;
//		}

	//--------------------- a template to get vectors from labels -------------------------
	//-------- Only use this once for now, potentially 13+ more lines can be saved --------


























//----------- Below are migrated from Singlephoton_module.cc
    void SinglePhoton::GetVertex(const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, const art::Ptr<recob::PFParticle> & particle ){

        if(m_is_verbose) std::cout<<"SinglePhoton::"<<__FUNCTION__<<"()\t||\t Starting to analyze recob::Vertex\n";
        int n_vert =0;

        //std::cout<<"There are "<<pfParticlesToVerticesMap.count(particle)<<" verticies associated with this particle"<<std::endl;

        lar_pandora::PFParticlesToVertices::const_iterator vIter = pfParticlesToVerticesMap.find(particle);
        if (pfParticlesToVerticesMap.end() != vIter)
        {
            const lar_pandora::VertexVector &vertexVector = vIter->second;
            if (!vertexVector.empty())
            {
                if (vertexVector.size() !=1)
                    std::cout << " Warning: Found particle with more than one associated vertex " << "\n";

                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);

                n_vert++;

                m_vertex_pos_x = xyz[0];
                m_vertex_pos_y = xyz[1];
                m_vertex_pos_z = xyz[2];
//				std::cout<<"CHECK Vertex position: ("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")"<<std::endl;
                std::vector<double> tmp = {xyz[0],xyz[1],xyz[2]};
                m_reco_vertex_in_SCB = this->distToSCB(m_reco_vertex_dist_to_SCB,tmp);
                m_reco_vertex_dist_to_active_TPC =  this->distToTPCActive(tmp);
//CHECK
//                if(!m_run_pi0_filter){
//                    m_reco_vertex_to_nearest_dead_wire_plane0 = distanceToNearestDeadWire(0, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
//                    m_reco_vertex_to_nearest_dead_wire_plane1 = distanceToNearestDeadWire(1, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
//                    m_reco_vertex_to_nearest_dead_wire_plane2 = distanceToNearestDeadWire(2, m_vertex_pos_y, m_vertex_pos_z,geom,bad_channel_list_fixed_mcc9);
//                }

            }else{
				throw art::Exception(art::errors::StdException)
                << " Pandor did not have an associated vertex for a particle. ";
//                << " Error: vertexVector associated with this particle is empty " << "\n";

            }
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Finished. Found "<<n_vert<<" vertices.\n";
    }

    void SinglePhoton::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap)
    {
        //std::cout<<"Filling pfParticleMap with from the handle with total number "<<pfParticleHandle->size()<<std::endl;
        for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
        {
            const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
            // std::cout<<"Adding PFP to pfParticleMap with pfp id  "<<pParticle->Self()<<std::endl;
            if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
            {
                throw cet::exception("SinglePhoton") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
            }
        }
    }



	//Classify PFParticles into crParticles or nuParticles.
	void SinglePhoton::GetFinalStatePFParticleVectors(
			const PFParticleIdMap &pfParticleMap, 
			const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, 
			PFParticleVector &crParticles, 
			PFParticleVector &nuParticles, 
			size_t fpfp_w_bestnuID ) {
		if(m_is_verbose) std::cout<<"SinglePhoton::"<<__FUNCTION__<<"\t||\tSort out PFPParticles."<<std::endl;

        int found = 0;
        int primaries = 0;
//CHECK		PFParticleIdMap::const_iterator best_nu = NULL;//only take the PFParticle with best nu score
        for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
        {
            const art::Ptr<recob::PFParticle> pParticle(it->second);

            // Only look for primary particles
            if (!pParticle->IsPrimary()) continue;
            primaries++;


            // Check if this particle is identified as the neutrino
            const int pdg(pParticle->PdgCode());
//            const bool isNeutrino =  (std::abs(pdg) ==  pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
			 //CHECK
             bool isNeutrino =  (std::abs(pdg) ==  12 || std::abs(pdg) == 14 || std::abs(pdg) == 16);


            // If it is, lets get the vertex position
            if(isNeutrino){
				if(found>0){
//					std::cout<<"CHECK WARNING we have more than 1 neutrinos here, not expected as in MicroBooNE but ok."<<std::endl;
				}
                found++;

//				std::cout<<"CHECK "<<pParticle->Self()<<"th PFParticle is reco. as neutrino w. pdg: "<<pdg<<std::endl;
				if(pParticle->Self() == fpfp_w_bestnuID){
					std::cout<<"Take this PFParticle's vertex as event vertex, bc of it has highest nu score"<<std::endl;
					this->GetVertex(pfParticlesToVerticesMap, pParticle );
				}
            }

            // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
            if (!isNeutrino)
            {
                crParticles.push_back(pParticle);
                continue;
            }

            // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
            //       If this is not the case please handle accordingly
//CHECK            if (!nuParticles.empty())
//            {
//                throw cet::exception("SinglePhoton") << "  This event contains "<<nuParticles.size()<<" reconstructed neutrinos! Not goot. ";
//            }

            // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
            for (const size_t daughterId : pParticle->Daughters())
            {
                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";

                nuParticles.push_back(pfParticleMap.at(daughterId));
            }
        }
        std::cout<<"SinglePhoton::"<<__FUNCTION__<<"\t||\t Found "<<primaries<<" primary PFParticles (out of "<<pfParticleMap.size()<<") of which: "<<found<<" were neutrinos.\n"<<std::endl;
        m_reco_vertex_size = found;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

	void SinglePhoton::CollectTracksAndShowers(
			const PFParticleVector &particles,
			const PFParticleIdMap pfParticleMap, 
			const PFParticleHandle &pfParticleHandle, 
			const art::Event &evt, 
			TrackVector &tracks, 
			ShowerVector &showers,  
			std::map< art::Ptr<recob::Track> , 
			art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, 
			std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
    {
        // Get the associations between PFParticles and tracks/showers from the event
        art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
        art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);

        //if running over the neutrino slice only 
        if (m_run_all_pfps == false){ 
            for (const art::Ptr<recob::PFParticle> &pParticle : particles) {
                const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
                const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));

                FillTracksAndShowers(associatedTracks, associatedShowers, pParticle,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);
            }
        } else{ //if running over all slices
            std::cout<<"SinglePhoton\t||\tThe total number of PFP's in the map is "<<pfParticleMap.size()<<std::endl;
            //            std::cout<<"The total number of PFP's in the vector is "<< particles.size()<<std::endl;
            for (auto pair : pfParticleMap){
                const art::Ptr<recob::PFParticle> &pParticle = pair.second;

                const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
                const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));

                FillTracksAndShowers(associatedTracks, associatedShowers, pParticle,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);

            }
        }
    }



    void SinglePhoton::FillTracksAndShowers( const std::vector< art::Ptr<recob::Track> > & associatedTracks, const std::vector< art::Ptr<recob::Shower> > & associatedShowers, const art::Ptr<recob::PFParticle> &pParticle , const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
    {

        const unsigned int nTracks(associatedTracks.size());
        const unsigned int nShowers(associatedShowers.size());


        // Check if the PFParticle has no associated tracks or showers
        if (nTracks == 0 && nShowers == 0)
        {
            //  std::cout<<"ERROR No tracks or showers were associated to PFParticle " << pParticle->Self()<<" with pdg "<<pParticle->PdgCode() <<std::endl;
            //std::cout<<"-- isPrimary = "<<pParticle->IsPrimary()<<std::endl;
            mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
            return;
        }

        // Check if there is an associated track
        if (nTracks == 1 && nShowers == 0)
        {

            tracks.push_back(associatedTracks.front());
            trackToNuPFParticleMap[tracks.back()]= pParticle;
            //std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  associatedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

            return;
        }

        // Check if there is an associated shower
        if (nTracks == 0 && nShowers == 1)
        {
            showers.push_back(associatedShowers.front());
            showerToNuPFParticleMap[showers.back()] = pParticle;
            // std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  associatedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

            return;
        }

        throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();

    }





    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input){
        corrected.resize(3);

//CHECK        double kx = input[0];
        double ky = input[1];
        double kz = input[2];
//CHECK
//        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
	   // CHECK
       // double g4Ticks =  detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->trigger_offset();

//CHECK       double xtimeoffset = 0;//CHECK  theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        //        double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = 0;//CHECK scecorr.Y();
        double zOffset = 0;//CHECK scecorr.Z();

        corrected[0]=0;//CHECK kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks->TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector->GetXTicksOffset(0,0,0)<<" "<<theDetector->trigger_offset()<<std::endl;
        return 0;
    }





    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);

        double kx = mcparticle->Vx();
        double ky = mcparticle->Vy();
        double kz = mcparticle->Vz();

//CHECK        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));  // to get position offsets to be used in ionization electron drift
//CHECK         double g4Ticks =  detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->trigger_offset();

//CHECK        double xtimeoffset =  theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        //double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = 0;//CHECK scecorr.Y();
        double zOffset = 0;//CHECK scecorr.Z();

        corrected[0]= kx;//CHECK- scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks->TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector->GetXTicksOffset(0,0,0)<<" "<<theDetector->trigger_offset()<<std::endl;
        return 0;
    }





    int SinglePhoton::spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);
        //Space Charge Effect! functionize this soon.
        double kx = mcparticle.Vx();
        double ky = mcparticle.Vy();
        double kz = mcparticle.Vz();
//CHECK         auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
//CHECK        double g4Ticks =  detClocks->TPCG4Time2Tick(mcparticle.T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->trigger_offset();

//CHECK        double xtimeoffset = 0;//CHECK theDetector->ConvertTicksToX(g4Ticks,0,0,0);

        corrected[0]=kx;//CHECK  - scecorr.X() +xtimeoffset+0.6;
        corrected[1]=ky;//CHECK  + scecorr.Y();
        corrected[2]=kz;//CHECK  + scecorr.Z();
        return 0;
    }

    void SinglePhoton::CollectMCParticles(
		const art::Event &evt,
		const std::string &label,
		std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,
		std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> &particlesToTruth,
		std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap){

        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

        art::Handle< std::vector< simb::MCParticle>  > theParticles;
        evt.getByLabel(label, theParticles);

        if (!theParticles.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
        }

        art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);

        for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
        {
            const art::Ptr<simb::MCParticle> particle(theParticles, i);
            const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
            truthToParticles[truth].push_back(particle);
            particlesToTruth[particle] = truth;
            MCParticleToTrackIdMap[particle->TrackId()] = particle;
        }

//        std::cout<<"SinglePhoton::CollectMCParticles() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
    }

    void SinglePhoton::CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector)
    {
        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectSimChannels --- Trying to access MC truth from real data ";

        art::Handle< std::vector<sim::SimChannel> > theSimChannels;
        evt.getByLabel(label, theSimChannels);

        if (!theSimChannels.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find sim channels... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theSimChannels->size() << " SimChannels " << std::endl;
        }

        for (unsigned int i = 0; i < theSimChannels->size(); ++i)
        {
            const art::Ptr<sim::SimChannel> channel(theSimChannels, i);
            simChannelVector.push_back(channel);
        }
    }


	void SinglePhoton::BuildMCParticleHitMaps(
			const art::Event &evt, 
			const std::string &label, 
			const std::vector<art::Ptr<recob::Hit>> &hitVector,   
			std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         
			std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, 
			const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, 
			std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap){
        std::vector< art::Ptr<sim::SimChannel> >   simChannelVector;
        std::map< art::Ptr<simb::MCTruth>,     std::vector<art::Ptr<simb::MCParticle>>  >    truthToParticles;
        std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> > particlesToTruth;
        std::map< art::Ptr<recob::Hit>,    std::vector< sim::TrackIDE >    >               hitsToTrackIDEs;

        this->CollectSimChannels(evt, label, simChannelVector);
        this->CollectMCParticles(evt, label, truthToParticles, particlesToTruth, MCParticleToTrackIdMap);
//        CHECK
		//Collect the links from reconstructed hits to their true energy deposits.
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(evt, hitVector, simChannelVector, hitsToTrackIDEs);
		//Build mapping between Hits and MCParticles, starting from Hit/TrackIDE/MCParticle information
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);


    }

    bool SinglePhoton::Pi0PreselectionFilter()
    {

        if(m_vertex_pos_x < m_tpc_active_XMin||  m_vertex_pos_x > m_tpc_active_XMax) return false;
        if(m_vertex_pos_y < m_tpc_active_YMin || m_vertex_pos_y > m_tpc_active_YMax) return false;
        if(m_vertex_pos_z < m_tpc_active_ZMin||  m_vertex_pos_z > m_tpc_active_ZMax) return false;

        if(m_reco_asso_showers!=2) return false;
        if(m_reco_asso_tracks!=1) return false;
        if(m_reco_vertex_size<1) return false;

        if(m_reco_shower_conversion_distance.size()!=2) return false;
        if(m_reco_shower_conversion_distance[0]<1. || m_reco_shower_conversion_distance[1]<1.) return false;

        return true;
    }



    bool SinglePhoton::Pi0PreselectionFilter2g0p()
    {
        if(m_vertex_pos_x < m_tpc_active_XMin||  m_vertex_pos_x > m_tpc_active_XMax) return false;
        if(m_vertex_pos_y < m_tpc_active_YMin || m_vertex_pos_y > m_tpc_active_YMax) return false;
        if(m_vertex_pos_z < m_tpc_active_ZMin||  m_vertex_pos_z > m_tpc_active_ZMax) return false;

        if(m_reco_asso_showers!=2) return false;
        if(m_reco_asso_tracks!=0) return false;
        if(m_reco_vertex_size<1) return false;

        if(m_reco_shower_energy_max.size()!=2) return false;
        //if the maximum energy of all showers on all planes is smaller than 30
        if(m_reco_shower_energy_max[m_reco_shower_ordered_energy_index[0]]<30.) return false;

        return true;
    }

    bool SinglePhoton::IsEventInList(int run, int subrun, int event){
        if(m_selected_set.find( {run, subrun, event} ) == m_selected_set.end()){
            if(m_selected_set.find({run, subrun})  == m_selected_set.end() ){
                if(m_selected_set.find({run}) == m_selected_set.end())
                    return false;
            }
        }
        return true;
    }

//----------- Above are migrated from Singlephoton_module.cc

    //determines if a point is inside the rectangle by summing the areas of the four triangles made by 
    //if the point is inside, the sum of the triangles should exactly equal the area of the rectangle
    //also returns true if the point is on the boundary
    bool SinglePhoton::isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle){
        int n_vertices = (int)rectangle.size();
        //bool inside = false;
        int i, j = 0;
        double areas = 0;
        //for each pair of vertices
        for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
            //calculate the area of a triangle with the point and two vertices
            double this_area = 0.5*abs(rectangle[i][0]*(rectangle[j][1] - thishit_pos[1]) 
									+ rectangle[j][0]*(thishit_pos[1] - rectangle[i][1]) 
									+ thishit_pos[0]*(rectangle[i][1] - rectangle[j][1]));
            areas += this_area;
        }        
        //calc area of the rectangle
        double area_rectangle = m_width_dqdx_box* m_length_dqdx_box;

        //check the sum of areas match
        if (abs(areas - area_rectangle) <= 0.001 ){
            return true;
        }
        return false;
    }


	//helpers for calculating calometry
    double SinglePhoton::CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits){
        double energy[3] = {0., 0., 0.};

        //std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t Looking at shower with "<<hits.size() <<" hits on all planes"<<std::endl;

        //for each hit in the shower
        for (auto &thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes     	
            if (plane > 2 || plane < 0)	continue;

            //calc the energy of the hit
            double E = QtoEConversion(GetQHit(thishitptr, plane));

            //add the energy to the plane
            energy[plane] += E;
        }//for each hiti

        //find the max energy on a single plane
        double max = energy[0];
        for (double en: energy){
            if( en > max){
                max = en;
            }
        }
        // std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The energy on each plane for this shower is "<<energy[0]<<", "<<energy[1]<<", "<<energy[2]<<std::endl;

        //return the highest energy on any of the planes
        return max;

    }

    double SinglePhoton::CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int this_plane){
        double energy = 0.;

        //for each hit in the shower
        for (auto &thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes     	
            if (plane != this_plane )	continue;

            //calc the energy of the hit
            double E = QtoEConversion(GetQHit(thishitptr, plane));
            //add the energy to the plane
            energy += E;
        }//for each hit

        return energy;

    }




    double SinglePhoton::GetQHit(art::Ptr<recob::Hit> thishitptr, int plane){
        double gain;
        //choose gain based on whether data/mc and by plane
        if (m_is_data == false &&  m_is_overlayed == false){
            gain = m_gain_mc[plane] ;
            //if (m_is_verbose) std::cout<<"the gain for mc on plane "<<plane<<" is "<<gain<<std::endl;
        } if (m_is_data == true ||  m_is_overlayed == true){
            gain = m_gain_data[plane] ;
            //if (m_is_verbose) std::cout<<"the gain for data on plane "<<plane<<" is "<<gain<<std::endl;

        }

        double Q = thishitptr->Integral()*gain;
        return Q;
    }


    double SinglePhoton::QtoEConversion(double Q){
        //return the energy value converted to MeV (the factor of 1e-6)
        double E = Q* m_work_function *1e-6 /m_recombination_factor;
        return E;

    }


    std::vector<double> SinglePhoton::CalcdEdxFromdQdx(std::vector<double> dqdx){
        int n = dqdx.size();
        std::vector<double> dedx(n,0.0);
        for (int i = 0; i < n; i++){
            //std::cout<<"The dQ/dx is "<<dqdx[i]<<std::endl;
            dedx[i] = QtoEConversion(dqdx[i]);
            //std::cout<<"The dE/dx is "<<dedx[i]<<std::endl;
        }
        return dedx;
    }


    std::vector<double> SinglePhoton::CalcdQdxShower(
            const art::Ptr<recob::Shower>& shower,
            const std::vector<art::Ptr<recob::Cluster>> & clusters, 
            std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane,
			double triggeroffset,
			detinfo::DetectorPropertiesData const & theDetector){
        //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The number of clusters in this shower is "<<clusters.size()<<std::endl;
        std::vector<double> dqdx;

        //get the 3D shower direction
        //note: in previous versions of the code there was a check to fix cases where the shower direction was inverted - this hasn't been implemented
        TVector3 shower_dir(shower->Direction().X(), shower->Direction().Y(),shower->Direction().Z());

        //calculate the pitch for this plane
        double pitch = getPitch(shower_dir, plane);	
        //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The pitch between the shower and plane "<<plane<<" is "<<pitch<<std::endl;

        //for all the clusters in the shower
        for (const art::Ptr<recob::Cluster> &thiscluster: clusters){
            //keep only clusters on the plane
            if(thiscluster->View() != plane) continue;

            //calculate the cluster direction
            std::vector<double> cluster_axis = {cos(thiscluster->StartAngle()), sin(thiscluster->StartAngle())};		

            //get the cluster start and and in CM
            //std::cout<<"for plane/tpc/cryo:"<<plane<<"/"<<m_TPC<<"/"<<m_Cryostat<<", fXTicksOffset: "<<theDetector->GetXTicksOffset(plane, m_TPC, m_Cryostat)<<" fXTicksCoefficient: "<<theDetector->GetXTicksCoefficient(m_TPC, m_Cryostat)<<std::endl;

            //convert the cluster start and end positions to time and wire coordinates
            std::vector<double> cluster_start = {thiscluster->StartWire() * m_wire_spacing,(thiscluster->StartTick() - triggeroffset)* _time2cm};
            std::vector<double> cluster_end = {thiscluster->EndWire() * m_wire_spacing,(thiscluster->EndTick() - triggeroffset)* _time2cm };

            //check that the cluster has non-zero length
            double length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) + pow(cluster_end[1] - cluster_start[1], 2));
            //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers() \t||\t The cluster length is "<<length<<std::endl;
            if (length <= 0){ 
                std::cout<<"skipping cluster on plane "<<plane<<", length = "<<length<<std::endl;
                continue;
            }


            //draw a rectangle around the cluster axis 
            std::vector<std::vector<double>> rectangle = buildRectangle(cluster_start, cluster_axis, m_width_dqdx_box, m_length_dqdx_box);	

            //get all the hits for this cluster
            std::vector<art::Ptr<recob::Hit>> hits =  clusterToHitMap[thiscluster];

            //for each hit in the cluster
            for (art::Ptr<recob::Hit> &thishit: hits){	
                //get the hit position in cm from the wire and time
                std::vector<double> thishit_pos ={thishit->WireID().Wire * m_wire_spacing, (thishit->PeakTime() - triggeroffset)* _time2cm};

                //check if inside the box
                bool v2 = isInsidev2(thishit_pos, rectangle);
                if (v2){
                    double q = GetQHit(thishit, plane); 
                    double this_dqdx = q/pitch; 
                    dqdx.push_back(this_dqdx);
                }//if hit falls inside the box

            }//for each hit inthe cluster
        }//for each cluster
        return dqdx;
    }

    double SinglePhoton::getPitch(TVector3 shower_dir, int plane){
        //get the wire direction for this plane - values are hardcoded which isn't great but the TPC geom object gave weird values
        TVector3 wire_dir = getWireVec(plane);

        //take dot product of shower and wire dir
        double cos = getCoswrtWires(shower_dir, wire_dir);

        //want only positive values so take abs, normalize by the lengths of the shower and wire
        cos = abs(cos)/(wire_dir.Mag() * shower_dir.Mag());	

        //If the cos is 0 shower is perpendicular and therefore get infinite distance 
        if (cos == 0){ return std::numeric_limits<double>::max(); }

        //output is always >= the wire spacing
        return m_wire_spacing/cos;
    }



    double SinglePhoton::getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
        int nhits = 0;
        double widths = 0;
        for (art::Ptr<recob::Hit> thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes       
            if (plane != this_plane) continue;

            widths += thishitptr->RMS(); // recob::Hit->RMS() returns RMS of the hit shape in tick units
            nhits++;
        }//for each hiti
        return   widths/(double)nhits;
    }



    int SinglePhoton::getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
        int nhits = 0;
        for (art::Ptr<recob::Hit> thishitptr : hits){
            //check the plane
            int plane= thishitptr->View();

            //skip invalid planes       
            if (plane != this_plane) continue;

            nhits++;
        }//for each hiti
        return nhits;
    }




    //area of a triangle given three vertices

//	double this_area = areaTriangle(rectangle[i][0], rectangle[i][1], rectangle[j][0], rectangle[j][1], thishit_pos[0], thishit_pos[1]);
//    double SinglePhoton::areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3){
//        double num = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);
//        return abs(num)/2;
//    }



//CHECK below are moved from SinglePhoton_module.cc


    double SinglePhoton::triangle_area(double a1, double a2, double b1, double b2, double c1, double c2){
        double m1 = 0.3;
        double m2 = 1.0/25.0;

        return fabs((a1*m1*(b2*m2-c2*m2)+b1*m1*(c2*m2-a2*m2)+c1*m1*(a2*m2-b2*m2))/2.0);
    }

    int SinglePhoton::quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area){

        std::vector<double> z(n,0.0);

        TGraph2D *g = new TGraph2D(n,X,Y,&z[0]);
        TGraphDelaunay delan(g);
        delan.SetMarginBinsContent(0);
        delan.ComputeZ(0,0);
        delan.FindAllTriangles();
        (*num_triangles)=delan.GetNdt(); // number of Delaunay triangles found

        //Grab the locations of all the trianges. These will be intergers referencing to position in X,Y arrays
        Int_t *MT = delan.GetMTried();
        Int_t *NT = delan.GetNTried();
        Int_t *PT = delan.GetPTried();

        (*area)=0.0;
        for(int i = 0; i<delan.GetNdt(); i++){
            (*area)+=triangle_area(X[MT[i]-1],Y[MT[i]-1],X[NT[i]-1],Y[NT[i]-1],X[PT[i]-1],Y[PT[i]-1]);
        }

        delete g;
        return 0;
    }

    int SinglePhoton::delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area){

        int n = hits.size();
        std::vector<double> C0,T0;
        std::vector<double> C1,T1;
        std::vector<double> C2,T2;
        size_t n_0=0;
        size_t n_1=0;
        size_t n_2=0;

        for(int i=0;i<n; i++){
            const art::Ptr<recob::Hit> hit = hits[i];
            switch(hit->View()){
                case 0:
                    C0.push_back((double)hit->WireID().Wire);         
                    T0.push_back(hit->PeakTime());         
                    n_0++;
                    break;
                case 1:
                    C1.push_back((double)hit->WireID().Wire);         
                    T1.push_back(hit->PeakTime());         
                    n_1++;
                    break;
                case 2:
                    C2.push_back((double)hit->WireID().Wire);         
                    T2.push_back(hit->PeakTime());         
                    n_2++;
                    break;
                default:
                    break;
            }
        }
        if(m_use_delaunay){
            if(n_0>0 && (int)n_0 < m_delaunay_max_hits) this->quick_delaunay_fit(n_0, &C0[0]  , &T0[0]  , &num_triangles[0],&area[0]);
            if(n_1>0 && (int)n_1 < m_delaunay_max_hits) this->quick_delaunay_fit(n_1, &C1[0]  , &T1[0]  , &num_triangles[1],&area[1]);
            if(n_2>0 && (int)n_2 < m_delaunay_max_hits) this->quick_delaunay_fit(n_2, &C2[0]  , &T2[0]  , &num_triangles[2],&area[2]);
        }
        num_hits[0] = n_0;
        num_hits[1] = n_1;
        num_hits[2] = n_2;

        //std::cout<<"Plane 0: "<<n_0<<" hits with "<<num_triangles[0]<<" triangles of area: "<< area[0]<<std::endl;
        //std::cout<<"Plane 1: "<<n_1<<" hits with "<<num_triangles[1]<<" triangles of area: "<< area[1]<<std::endl;
        //std::cout<<"Plane 2: "<<n_2<<" hits with "<<num_triangles[2]<<" triangles of area: "<< area[2]<<std::endl;

        return 0;
    }

}
