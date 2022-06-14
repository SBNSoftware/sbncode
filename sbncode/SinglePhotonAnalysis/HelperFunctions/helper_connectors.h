namespace single_photon
{
	//--------------------- a template to get vectors from labels -------------------------
	//1. Scratch from Handle, and return a equivalently useful vector.
	//  Sample usage:
	//		art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); //yea, this should be gone;
	//		recob::Hit dummy_hit;//This is to specify the template;
	//		std::vector<art::Ptr<recob::Hit>> hitVector = VectorFromLabel(dummy_hit, evt, m_hitfinderLabel);

	template <typename recob_object>//A helper template that allows you to make compliated types.
		struct temporary_types{ 
			using type1 = std::vector<art::Ptr<recob_object>>;
			using type2 = art::ValidHandle<std::vector<recob_object>>;
			using type3 = std::vector<recob_object>;
		};

	template <class recob_object>//ref_type is only used to identify the temporary_types late.
		typename temporary_types<recob_object>::type1 VectorFromLabel(recob_object ref_type, const art::Event &evt, std::string &label){

			typename temporary_types<recob_object>::type2 const & Handle = evt.getValidHandle<typename temporary_types<recob_object>::type3>(label);
			typename temporary_types<recob_object>::type1 Vector;
			art::fill_ptr_vector(Vector,Handle);
			return Vector;
		}

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
				std::cout<<"CHECK Vertex position: ("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")"<<std::endl;
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
    void SinglePhoton::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, PFParticleVector &crParticles, PFParticleVector &nuParticles, size_t fpfp_w_bestnuID )
    {
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
            //const bool isNeutrino =  (std::abs(pdg) ==   14 );


            // If it is, lets get the vertex position
            if(isNeutrino){
				if(found>0){
					std::cout<<"CHECK WARNING we have more than 1 neutrinos here, not expected as in MicroBooNE but ok."<<std::endl;
				}
                found++,

				std::cout<<"CHECK "<<pParticle->Self()<<"th PFParticle is reco. as neutrino w. pdg: "<<pdg<<std::endl;
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

    void SinglePhoton::CollectTracksAndShowers(const PFParticleVector &particles,const PFParticleIdMap pfParticleMap, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
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

        std::cout<<"SinglePhoton::CollectMCParticles() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
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

}
