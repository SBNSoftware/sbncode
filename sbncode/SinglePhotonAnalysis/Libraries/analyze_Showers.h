namespace single_photon
{
	void SinglePhoton::AnalyzeShowers(
			std::vector<PandoraPFParticle> all_PPFPs,
			const std::vector<art::Ptr<recob::Shower>>& showers,  
//			std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap, 
//			std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap, 
//			std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > & pfParticleToClusterMap,
			std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> >  & clusterToHitMap , 
//			std::map<int, double>& sliceIdToNuScoreMap,
//			std::map<art::Ptr<recob::PFParticle>,bool>& PFPToClearCosmicMap,
//			std::map<art::Ptr<recob::PFParticle>, int>& PFPToSliceIdMap, 
//			std::map<art::Ptr<recob::PFParticle>,bool> &PFPToNuSliceMap, 
//			std::map<art::Ptr<recob::PFParticle>,double> &PFPToTrackScoreMap,
//			PFParticleIdMap &pfParticleMap,
			std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>>& PFPtoShowerReco3DMap,
			double triggeroffset,
			detinfo::DetectorPropertiesData const & theDetector
			){
//CHECK try to use all_PPFPs to get the correct nu score.
		m_is_verbose= true;
//        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Begininning recob::Shower analysis suite"<<std::endl;;

        m_reco_asso_showers=showers.size();
        int i_shr = 0;
        this->ResizeShowers(m_reco_asso_showers);
		
		if(m_is_verbose) std::cout<<__FUNCTION__<<" CHECK showers found: "<<showers.size()<<std::endl;
//		std::vector<int> spacers = Printer_header({"Slice","pfpID","Start(x,  ","   y,      ",",      z  )"});
        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
        {


            const art::Ptr<recob::Shower> shower = *iter;
//            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];
			//CHECK I think I can take care of all pfpmaps here;
			PandoraPFParticle* ppfp = PPFP_GetPPFPFromShower(all_PPFPs, shower);

			const art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;
//			std::cout<<"CHECK SHower pfparticle "<<pfp->Self()<<" and "<<pfp2->Self()<<std::endl;

            art::Ptr<recob::Shower> shower3d;
            if(PFPtoShowerReco3DMap.count(pfp)==0){
//                std::cout<<"SHOWER_3D : <<CHECK!! No shower3d in map for this pfp"<<std::endl;
//                std::cout<<"Reverting to normal recob::Shower"<<std::endl;
                m_reco_shower3d_exists[i_shr] = 0;
                shower3d = shower;
            }else{
                shower3d   = PFPtoShowerReco3DMap[pfp];
                m_reco_shower3d_exists[i_shr] = 1;
            }
//            const std::vector<art::Ptr<recob::Hit>> hits =  pfParticleToHitMap[pfp];
//			const std::vector<art::Ptr<recob::Cluster>> clusters = pfParticleToClusterMap[pfp];

            const std::vector<art::Ptr<recob::Hit>> hits =  ppfp->pHits;
			const std::vector<art::Ptr<recob::Cluster>> clusters = ppfp->pClusters;

            //int m_shrid = shower->ID(); This is an used variable, always -999
            double m_length = shower->Length();
            double m_open_angle = shower->OpenAngle();

            TVector3 shr_start = shower->ShowerStart();
            TVector3 shr_dir = shower->Direction();

            TVector3 shr3d_start = shower3d->ShowerStart();
            TVector3 shr3d_dir = shower3d->Direction();

//            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t On Shower: "<<i_shr<<" which has length: "<<m_length<<""<<std::endl;;

            m_reco_shower_startx[i_shr] = shr_start.X();
            m_reco_shower_starty[i_shr] = shr_start.Y();
            m_reco_shower_startz[i_shr] = shr_start.Z();

 
            std::vector<double> hstart = {m_reco_shower_startx[i_shr],m_reco_shower_starty[i_shr],m_reco_shower_startz[i_shr]};
            m_reco_shower_start_dist_to_active_TPC[i_shr] = distToTPCActive(hstart);
            m_reco_shower_start_in_SCB[i_shr] = this->distToSCB(m_reco_shower_start_dist_to_SCB[i_shr],hstart);

            m_reco_shower_dirx[i_shr] = shr_dir.X();
            m_reco_shower_diry[i_shr] = shr_dir.Y();
            m_reco_shower_dirz[i_shr] = shr_dir.Z();
            m_reco_shower_length[i_shr] = m_length;
            m_reco_shower_openingangle[i_shr] = m_open_angle;

            m_reco_shower3d_startx[i_shr] = shr3d_start.X();
            m_reco_shower3d_starty[i_shr] = shr3d_start.Y();
            m_reco_shower3d_startz[i_shr] = shr3d_start.Z();
            m_reco_shower3d_dirx[i_shr] = shr3d_dir.X();
            m_reco_shower3d_diry[i_shr] = shr3d_dir.Y();
            m_reco_shower3d_dirz[i_shr] = shr3d_dir.Z();
            m_reco_shower3d_length[i_shr] = shower3d->Length();
            m_reco_shower3d_openingangle[i_shr] = shower3d->OpenAngle();


            m_reco_shower_conversion_distance[i_shr] = sqrt( pow(shr_start.X()-m_vertex_pos_x,2)+pow(shr_start.Y()-m_vertex_pos_y,2)+ pow(shr_start.Z()-m_vertex_pos_z,2)  );
            m_reco_shower3d_conversion_distance[i_shr] = sqrt( pow(shr3d_start.X()-m_vertex_pos_x,2)+pow(shr3d_start.Y()-m_vertex_pos_y,2)+ pow(shr3d_start.Z()-m_vertex_pos_z,2)  );

            //pandroa shower
            std::vector<double> shr_ts = {shr_start.X(), shr_start.Y(), shr_start.Z()};
            std::vector<double> shr_te = {shr_start.X()-shr_dir.X(),shr_start.Y()-shr_dir.Y(),shr_start.Z()-shr_dir.Z()};
            std::vector<double> shr_tv = {m_vertex_pos_x,m_vertex_pos_y,m_vertex_pos_z};

            m_reco_shower_impact_parameter[i_shr] = dist_line_point(shr_ts,shr_te,shr_tv );
            m_reco_shower_implied_dirx[i_shr] = shr_start.X()-m_vertex_pos_x;;
            m_reco_shower_implied_diry[i_shr] = shr_start.Y()-m_vertex_pos_y;
            m_reco_shower_implied_dirz[i_shr] = shr_start.Z()-m_vertex_pos_z;

            double norm = sqrt(pow(m_reco_shower_implied_dirx[i_shr],2)+pow(m_reco_shower_implied_diry[i_shr],2)+pow(m_reco_shower_implied_dirz[i_shr],2));
            m_reco_shower_implied_dirx[i_shr] = m_reco_shower_implied_dirx[i_shr]/norm;
            m_reco_shower_implied_diry[i_shr] = m_reco_shower_implied_diry[i_shr]/norm;
            m_reco_shower_implied_dirz[i_shr] = m_reco_shower_implied_dirz[i_shr]/norm;

            //now 3D shower
            std::vector<double> shr3d_ts = {shr3d_start.X(), shr3d_start.Y(), shr3d_start.Z()};
            std::vector<double> shr3d_te = {shr3d_start.X()-shr3d_dir.X(),shr3d_start.Y()-shr3d_dir.Y(),shr3d_start.Z()-shr3d_dir.Z()};
            std::vector<double> shr3d_tv = {m_vertex_pos_x,m_vertex_pos_y,m_vertex_pos_z};

            m_reco_shower3d_impact_parameter[i_shr] = dist_line_point(shr3d_ts,shr3d_te,shr3d_tv );
            m_reco_shower3d_implied_dirx[i_shr] = shr3d_start.X()-m_vertex_pos_x;;
            m_reco_shower3d_implied_diry[i_shr] = shr3d_start.Y()-m_vertex_pos_y;
            m_reco_shower3d_implied_dirz[i_shr] = shr3d_start.Z()-m_vertex_pos_z;

            double shr3d_norm = sqrt(pow(m_reco_shower3d_implied_dirx[i_shr],2)+pow(m_reco_shower3d_implied_diry[i_shr],2)+pow(m_reco_shower3d_implied_dirz[i_shr],2));
            m_reco_shower3d_implied_dirx[i_shr] = m_reco_shower3d_implied_dirx[i_shr]/shr3d_norm;
            m_reco_shower3d_implied_diry[i_shr] = m_reco_shower3d_implied_diry[i_shr]/shr3d_norm;
            m_reco_shower3d_implied_dirz[i_shr] = m_reco_shower3d_implied_dirz[i_shr]/shr3d_norm;


            m_reco_shower_theta_yz[i_shr] = atan2(m_reco_shower_diry[i_shr],m_reco_shower_dirz[i_shr]);
            m_reco_shower_phi_yx[i_shr] = atan2(m_reco_shower_diry[i_shr],m_reco_shower_dirx[i_shr]);

            m_reco_shower3d_theta_yz[i_shr] = atan2(m_reco_shower3d_diry[i_shr],m_reco_shower3d_dirz[i_shr]);
            m_reco_shower3d_phi_yx[i_shr] = atan2(m_reco_shower3d_diry[i_shr],m_reco_shower3d_dirx[i_shr]);


            m_reco_shower_start_to_nearest_dead_wire_plane0[i_shr] = distanceToNearestDeadWire(0, m_reco_shower_starty[i_shr], m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
            m_reco_shower_start_to_nearest_dead_wire_plane1[i_shr] = distanceToNearestDeadWire(1, m_reco_shower_starty[i_shr], m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
            m_reco_shower_start_to_nearest_dead_wire_plane2[i_shr] = distanceToNearestDeadWire(2, m_reco_shower_starty[i_shr], m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
            std::vector<int> t_num(3,0);   // num of triangles on each plane
            std::vector<int> t_numhits(3,0);  // num of hits on each plane
            std::vector<double> t_area(3,0.0);

            //Right, this basically loops over all hits in all planes and for each plane forms the Delaunay triangilization of it and calculates the 2D area inscribed by the convex hull
//            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Starting Delaunay Triangleization"<<std::endl;;

            //auto start = std::chrono::high_resolution_clock::now();
            this->delaunay_hit_wrapper(hits, t_numhits, t_num, t_area);

            //auto finish = std::chrono::high_resolution_clock::now();
            //auto microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
            //if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t Finished Delaunay Triangleization. It took "<< microseconds.count() << "ms and found "<<t_num[0]+t_num[1]+t_num[2]<<" triangles"<<std::endl;;

            m_reco_shower_delaunay_num_triangles_plane0[i_shr] = t_num[0];
            m_reco_shower_delaunay_num_triangles_plane1[i_shr] = t_num[1];
            m_reco_shower_delaunay_num_triangles_plane2[i_shr] = t_num[2];

            m_reco_shower_delaunay_area_plane0[i_shr] = t_area[0];
            m_reco_shower_delaunay_area_plane1[i_shr] = t_area[1];
            m_reco_shower_delaunay_area_plane2[i_shr] = t_area[2];

            m_reco_shower_num_hits_plane0[i_shr] = t_numhits[0];
            m_reco_shower_num_hits_plane1[i_shr] = t_numhits[1];
            m_reco_shower_num_hits_plane2[i_shr] = t_numhits[2];
            //-------------- Calorimetry 3D --------------------


            const std::vector< double > shr3d_energy = shower3d->Energy();
            const std::vector< double > shr3d_dEdx = shower3d->dEdx();
            //const int shr3d_bestplane = shower3d->best_plane();

 //           std::cout<<"SHOWER3D_ENERGY: best plane: "<<shr3d_bestplane<<std::endl;
            //for(auto &en:shr3d_energy){
            //    std::cout<<en<<" ";
            //}
            if(shr3d_energy.size()==3){
                m_reco_shower3d_energy_plane0[i_shr] = shr3d_energy[0];
                m_reco_shower3d_energy_plane1[i_shr] = shr3d_energy[1];
                m_reco_shower3d_energy_plane2[i_shr] = shr3d_energy[2];
            }else{
                m_reco_shower3d_energy_plane0[i_shr] =-99;
                m_reco_shower3d_energy_plane1[i_shr] =-99;
                m_reco_shower3d_energy_plane2[i_shr] =-999;
            }

   //         std::cout<<std::endl<<"SHOWER3D_DEDX: "<<std::endl;
            //for(auto &dedx: shr3d_dEdx){
            //    std::cout<<dedx<<" ";
            //}
            if(shr3d_dEdx.size()==3){
                m_reco_shower3d_dEdx_plane0[i_shr] = shr3d_dEdx[0];
                m_reco_shower3d_dEdx_plane1[i_shr] = shr3d_dEdx[1];
                m_reco_shower3d_dEdx_plane2[i_shr] = shr3d_dEdx[2];
            }else{
                m_reco_shower3d_dEdx_plane0[i_shr] =-99;
                m_reco_shower3d_dEdx_plane1[i_shr] =-99;
                m_reco_shower3d_dEdx_plane2[i_shr] =-999;
            }


            //------------- calorimetry ------------

            m_reco_shower_energy_max[i_shr] = CalcEShower(hits);
            m_reco_shower_energy_plane0[i_shr] = CalcEShowerPlane(hits, 0);
            m_reco_shower_energy_plane1[i_shr] = CalcEShowerPlane(hits, 1);
            m_reco_shower_energy_plane2[i_shr] = CalcEShowerPlane(hits, 2);

            m_reco_shower_plane0_nhits[i_shr] = getNHitsPlane(hits, 0);
            m_reco_shower_plane1_nhits[i_shr] = getNHitsPlane(hits, 1);
            m_reco_shower_plane2_nhits[i_shr] = getNHitsPlane(hits, 2);

            m_reco_shower_plane0_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 0);
            m_reco_shower_plane1_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 1);
            m_reco_shower_plane2_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 2);


            //currently only run on 1 shower events
            if(showers.size()==1){
               for(auto &h: hits){ 

                    int plane= h->View();
                    int wire = h->WireID().Wire;
                    int tick = h->PeakTime();

                    m_reco_shower_hit_tick.push_back(tick);
                    m_reco_shower_hit_plane.push_back(plane);
                    m_reco_shower_hit_wire.push_back(wire);




               }

            }



            //std::cout<<"The energy on each plane is 0: "<< m_reco_shower_energy_plane0[i_shr]<<", 1: "<< m_reco_shower_energy_plane1[i_shr]<<", 2: "<<  m_reco_shower_energy_plane2[i_shr]<<std::endl;


//			if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 0"<<std::endl;
//			if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 1"<<std::endl;
//			if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t starting dq/dx plane 2"<<std::endl; 
//
//			if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tstarting de/dx plane 0"<<std::endl;
//			if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tstarting de/dx plane 1"<<std::endl;
//			if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tstarting de/dx plane 2"<<std::endl;
//            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tPAR: "<<m_reco_shower_dEdx_plane0[i_shr].size()<<" "<<m_reco_shower_dEdx_plane1[i_shr].size()<<" "<<m_reco_shower_dEdx_plane2[i_shr].size()<<std::endl;

			m_reco_shower_dQdx_plane0[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 0 , triggeroffset, theDetector);
            m_reco_shower_dQdx_plane1[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 1 , triggeroffset, theDetector);
            m_reco_shower_dQdx_plane2[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 2 , triggeroffset, theDetector);
            m_reco_shower_dEdx_plane0[i_shr] = CalcdEdxFromdQdx(m_reco_shower_dQdx_plane0[i_shr]);
            m_reco_shower_dEdx_plane1[i_shr] = CalcdEdxFromdQdx(m_reco_shower_dQdx_plane1[i_shr]);

            m_reco_shower_dEdx_plane2[i_shr] = CalcdEdxFromdQdx(m_reco_shower_dQdx_plane2[i_shr]);

            m_reco_shower_dEdx_plane0_median[i_shr] = getMedian(m_reco_shower_dEdx_plane0[i_shr]);
            m_reco_shower_dEdx_plane1_median[i_shr] = getMedian(m_reco_shower_dEdx_plane1[i_shr]);
            m_reco_shower_dEdx_plane2_median[i_shr] = getMedian(m_reco_shower_dEdx_plane2[i_shr]);

            m_reco_shower_angle_wrt_wires_plane0[i_shr] = getAnglewrtWires(shr_dir,0);
            m_reco_shower_angle_wrt_wires_plane1[i_shr] = getAnglewrtWires(shr_dir,1);
            m_reco_shower_angle_wrt_wires_plane2[i_shr] = getAnglewrtWires(shr_dir,2);


            m_reco_shower_dQdx_plane0_median[i_shr] = getMedian(m_reco_shower_dQdx_plane0[i_shr]);
            m_reco_shower_dQdx_plane1_median[i_shr] = getMedian(m_reco_shower_dQdx_plane1[i_shr]);
            m_reco_shower_dQdx_plane2_median[i_shr] = getMedian(m_reco_shower_dQdx_plane2[i_shr]);



            m_reco_shower_dEdx_plane0_mean[i_shr] = std::accumulate(m_reco_shower_dEdx_plane0[i_shr].begin(), m_reco_shower_dEdx_plane0[i_shr].end(), 0.0)/((double)m_reco_shower_dEdx_plane0[i_shr].size()); 
            m_reco_shower_dEdx_plane1_mean[i_shr] = std::accumulate(m_reco_shower_dEdx_plane1[i_shr].begin(), m_reco_shower_dEdx_plane1[i_shr].end(), 0.0)/((double)m_reco_shower_dEdx_plane1[i_shr].size()); 
            m_reco_shower_dEdx_plane2_mean[i_shr] = std::accumulate(m_reco_shower_dEdx_plane2[i_shr].begin(), m_reco_shower_dEdx_plane2[i_shr].end(), 0.0)/((double)m_reco_shower_dEdx_plane2[i_shr].size()); 

            auto maxp0 = std::max_element(m_reco_shower_dEdx_plane0[i_shr].begin(), m_reco_shower_dEdx_plane0[i_shr].end());
            auto maxp1 = std::max_element(m_reco_shower_dEdx_plane1[i_shr].begin(), m_reco_shower_dEdx_plane1[i_shr].end());
            auto maxp2 = std::max_element(m_reco_shower_dEdx_plane2[i_shr].begin(), m_reco_shower_dEdx_plane2[i_shr].end());
            auto minp0 = std::min_element(m_reco_shower_dEdx_plane0[i_shr].begin(), m_reco_shower_dEdx_plane0[i_shr].end());
            auto minp1 = std::min_element(m_reco_shower_dEdx_plane1[i_shr].begin(), m_reco_shower_dEdx_plane1[i_shr].end());
            auto minp2 = std::min_element(m_reco_shower_dEdx_plane2[i_shr].begin(), m_reco_shower_dEdx_plane2[i_shr].end());


            if(maxp0 == m_reco_shower_dEdx_plane0[i_shr].end()){
                m_reco_shower_dEdx_plane0_max[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane0_max[i_shr] = *maxp0; 
            }

            if(maxp1 == m_reco_shower_dEdx_plane1[i_shr].end()){
                m_reco_shower_dEdx_plane1_max[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane1_max[i_shr] = *maxp1; 
            }

            if(maxp2 == m_reco_shower_dEdx_plane2[i_shr].end()){
                m_reco_shower_dEdx_plane2_max[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane2_max[i_shr] = *maxp2; 
            }


            if(minp0 == m_reco_shower_dEdx_plane0[i_shr].end()){
                m_reco_shower_dEdx_plane0_min[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane0_min[i_shr] = *minp0; 
            }

            if(minp1 == m_reco_shower_dEdx_plane1[i_shr].end()){
                m_reco_shower_dEdx_plane1_min[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane1_min[i_shr] = *minp1; 
            }

            if(minp2 == m_reco_shower_dEdx_plane2[i_shr].end()){
                m_reco_shower_dEdx_plane2_min[i_shr] = -999; 
            }else{
                m_reco_shower_dEdx_plane2_min[i_shr] = *minp2; 
            }


            m_reco_shower_dEdx_plane0_nhits[i_shr] = m_reco_shower_dEdx_plane0[i_shr].size();
            m_reco_shower_dEdx_plane1_nhits[i_shr] = m_reco_shower_dEdx_plane1[i_shr].size();
            m_reco_shower_dEdx_plane2_nhits[i_shr] = m_reco_shower_dEdx_plane2[i_shr].size();

            m_reco_shower_dEdx_amalgamated[i_shr] = getAmalgamateddEdx( 
					m_reco_shower_angle_wrt_wires_plane0[i_shr],  
					m_reco_shower_angle_wrt_wires_plane1[i_shr],  
					m_reco_shower_angle_wrt_wires_plane2[i_shr], 
					m_reco_shower_dEdx_plane0_median[i_shr], 
					m_reco_shower_dEdx_plane1_median[i_shr], 
					m_reco_shower_dEdx_plane2_median[i_shr],
					m_reco_shower_dEdx_plane0_nhits[i_shr], 
					m_reco_shower_dEdx_plane1_nhits[i_shr], 
					m_reco_shower_dEdx_plane2_nhits[i_shr] );

			m_reco_shower_dEdx_amalgamated_nhits[i_shr] = getAmalgamateddEdxNHits(
					m_reco_shower_dEdx_amalgamated[i_shr], 
					m_reco_shower_dEdx_plane0_median[i_shr], 
					m_reco_shower_dEdx_plane1_median[i_shr], 
					m_reco_shower_dEdx_plane2_median[i_shr],
					m_reco_shower_dEdx_plane0_nhits[i_shr], 
					m_reco_shower_dEdx_plane1_nhits[i_shr], 
					m_reco_shower_dEdx_plane2_nhits[i_shr] );

            //-------------- Flashes : Was there a flash in the beam_time and if so was it near in Z? --------------------
            double zmin = m_reco_shower_startz[i_shr];
            double zmax = zmin + m_reco_shower_dirz[i_shr]*m_reco_shower_length[i_shr];
            if(zmin > zmax) std::swap(zmin, zmax);

            double ymin = m_reco_shower_starty[i_shr];
            double ymax = zmin + m_reco_shower_diry[i_shr]*m_reco_shower_length[i_shr];
            if(ymin > ymax) std::swap(ymin, ymax);

            //Code property of Gray Yarbrough (all rights reserved)
            //int optical_flash_in_beamgate_counter=0;
            double shortest_dist_to_flash_z=DBL_MAX;
            double shortest_dist_to_flash_y=DBL_MAX;
            double shortest_dist_to_flash_yz=DBL_MAX;
            //-999 my nonsenese int can change
            int shortest_dist_to_flash_index_z=-999;
            int shortest_dist_to_flash_index_y=-999;
            int shortest_dist_to_flash_index_yz=-999;

//            if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\tnumber of flashes: "<< m_reco_num_flashes<< ""<<std::endl;;
            for(int i_flash = 0; i_flash < m_reco_num_flashes; ++i_flash) {

                double const zcenter=m_reco_flash_zcenter[i_flash];
//                if(m_is_verbose) std::cout<< "SinglePhoton::AnalyzeShowers()\t||\tflash z center:" <<m_reco_flash_zcenter[i_flash]<< ""<<std::endl;;
                double const ycenter=m_reco_flash_ycenter[i_flash];
//                if(m_is_verbose) std::cout<< "SinglePhoton::AnaluzeShowers()\t||\tflash y center:" <<m_reco_flash_ycenter[i_flash]<< ""<<std::endl;;

                //z plane
                double dist_z=DBL_MAX;
                if(zcenter < zmin) {
                    dist_z = zmin - zcenter;
                }
                else if(zcenter > zmax) {
                    dist_z = zcenter - zmax;
                }
                else {
                    dist_z = 0;
                }	    
                if(dist_z < shortest_dist_to_flash_z){
                    shortest_dist_to_flash_z = dist_z;
                    shortest_dist_to_flash_index_z=i_flash;
                }


                //y plane

                double dist_y=DBL_MAX;
                if(ycenter < ymin) {
                    dist_y = ymin - ycenter;
                }
                else if(ycenter > ymax) {
                    dist_y = ycenter - ymax;
                }
                else {
                    dist_y= 0;
                }	    
                if(dist_y < shortest_dist_to_flash_y){
                    shortest_dist_to_flash_y = dist_y;
                    shortest_dist_to_flash_index_y=i_flash;
                }

                double dist_yz=DBL_MAX;
                dist_yz=std::sqrt(dist_y*dist_y+dist_z*dist_z);
                if(dist_yz<shortest_dist_to_flash_yz){
                    shortest_dist_to_flash_yz = dist_yz;
                    shortest_dist_to_flash_index_yz=i_flash;
                }

            }


            //assume setting to nonsense value
            if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_z = -2;
            m_reco_shower_flash_shortest_distz[i_shr]=shortest_dist_to_flash_z;
            m_reco_shower_flash_shortest_index_z[i_shr]=shortest_dist_to_flash_index_z;

            if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_y = -2;
            m_reco_shower_flash_shortest_disty[i_shr]=shortest_dist_to_flash_y;
            m_reco_shower_flash_shortest_index_y[i_shr]=shortest_dist_to_flash_index_y;
            m_reco_shower_flash_shortest_distyz[i_shr]=shortest_dist_to_flash_yz;
            m_reco_shower_flash_shortest_index_yz[i_shr]=shortest_dist_to_flash_index_yz;
            if(m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_yz = -2;

//            if(m_is_verbose) std::cout << "the shortest distance in z plane between a flash and the shower: " << shortest_dist_to_flash_z << ""<<std::endl;;
//            if(m_is_verbose) std::cout << "the index closest flash to shower z_plane: " << shortest_dist_to_flash_index_z << ""<<std::endl;;
//            if(m_is_verbose) std::cout <<"the shortest distance in y plane between a flash and the shower: " << shortest_dist_to_flash_y << ""<<std::endl;;
//            if(m_is_verbose) std::cout << "the index closest flash to shower y_plane: " << shortest_dist_to_flash_index_y << ""<<std::endl;;
//            if(m_is_verbose) std::cout <<"the shortest distance in yz between a flash and the shower: " << shortest_dist_to_flash_yz << ""<<std::endl;;
//            if(m_is_verbose) std::cout << "the index closest flash to shower yz: " << shortest_dist_to_flash_index_yz << ""<<std::endl;;


            //end optical flash code


            m_reco_shower_num_daughters[i_shr] = pfp->NumDaughters();  //corresponding PFParticle
//			std::cout<<" CHECK numebr "<<m_reco_shower_num_daughters[i_shr]<<std::endl;
            if(m_reco_shower_num_daughters[i_shr]>0){
                //currently just look at 1 daughter
                //m_reco_shower_daughter_trackscore[i_shr] = PFPToTrackScoreMap[pfParticleMap[pfp->Daughters().front()]];
				int pfp_size = all_PPFPs.size();
				for(int index = 0; index < pfp_size; index++){
//					std::cout<<"CHECK Compare "<<pfp->Daughters().front()<<
//					" "<<all_PPFPs[index].pPFParticle->Self()<<std::endl;
					if( (pfp->Daughters().front()) == all_PPFPs[index].pPFParticle->Self());
					m_reco_shower_daughter_trackscore[i_shr] = all_PPFPs[index].get_TrackScore();
					break;
				}
            }


            //------------and finally some slice info-----------------

            m_reco_shower_sliceId[i_shr] = ppfp->get_SliceID();//PFPToSliceIdMap[pfp];
            m_reco_shower_nuscore[i_shr] = ppfp->get_NuScore();//sliceIdToNuScoreMap[ m_reco_shower_sliceId[i_shr]] ;
            m_reco_shower_isclearcosmic[i_shr] = ppfp->get_IsClearCosmic();//PFPToClearCosmicMap[pfp];
            m_reco_shower_is_nuslice[i_shr] = ppfp->get_IsNuSlice();//PFPToNuSliceMap[pfp];
			//m_reco_shower_trackscore[i_shr] = PFPToTrackScoreMap[pfp];
			//std::cout<<"m_reco_shower_is_nuslice[i_shr] = "<<m_reco_shower_is_nuslice[i_shr]<<" for shr with pfp "<<pfp->Self()<<std::endl; 


			m_reco_shower_trackscore[i_shr] = ppfp->get_TrackScore();
			m_reco_shower_pfparticle_pdg[i_shr] = ppfp->get_PdgCode();
           // if ( PFPToTrackScoreMap.find(pfp) != PFPToTrackScoreMap.end() ) {
           //     m_reco_shower_trackscore[i_shr] = PFPToTrackScoreMap[pfp];
           //     m_reco_shower_pfparticle_pdg[i_shr] = pfp->PdgCode();
           // } else{
           //     m_reco_shower_trackscore[i_shr] = -999; 
           //     m_reco_shower_pfparticle_pdg[i_shr] = -999;
           // }

//            if ( m_reco_shower_sliceId[i_shr] >0) std::cout<<"SinglePhoton::AnalyzeShowers()\t||\t On Shower: "<<i_shr<<". Pfp id = "<< pfp->Self()<<". The slice id for this shower is "<< m_reco_shower_sliceId[i_shr]<<", the neutrino score for this slice is "<< m_reco_shower_nuscore[i_shr]<<", and is_nuslice = "<<  m_reco_shower_is_nuslice[i_shr]<<". The track score is : "<< m_reco_shower_trackscore[i_shr]<<std::endl;

            i_shr++;

		//std::vector<int> spacers = Printer_header({"Slice","pfpID","Start(x,  ","   y,      ",",      z  )"});
//			Printer_content(
//					{std::to_string(ppfp->pSliceID),
//					std::to_string(ppfp->pPFParticleID),
//					std::to_string(shr_start.X()),
//					std::to_string(shr_start.Y()),
//					std::to_string(shr_start.Z())
//					},spacers);

        }

        //Lets sort and order the showers
        m_reco_shower_ordered_energy_index = sort_indexes(m_reco_shower_energy_max);
    }


	void SinglePhoton::AnalyzeKalmanShowers(
			const std::vector<art::Ptr<recob::Shower>>& showers,
			std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> &showerToPFParticleMap,
			std::map<art::Ptr<recob::PFParticle>,art::Ptr<recob::Track>> & pfParticlesToShowerKalmanMap,
			std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>>&  kalmanTrackToCaloMap, 
			std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap,
			detinfo::DetectorPropertiesData const & theDetector){
		std::cout<<"NO Kalman stuff in SBND" <<std::endl;
		return;
//KENG
//        std::cout<<"Singlephoton::"<<__FUNCTION__<<"\t||\tStarting to Analyze  (n="<<showers.size()<<")a showers via Kalman "<<std::endl;
//
//        std::vector<double> gains = {0,0,0};
//        for(int plane =0; plane < 3; plane++){
//            if (m_is_data == false &&  m_is_overlayed == false){
//                gains[plane] = m_gain_mc[plane] ;
//            } if (m_is_data == true ||  m_is_overlayed == true){
//                gains[plane] = m_gain_data[plane] ;
//            }
//        }
//
//
//        int i_shr=0;
//        for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
//        {
//            const art::Ptr<recob::Shower> shower = *iter;
//            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];
//            std::vector<art::Ptr<recob::Hit>> hitz = pfParticleToHitMap[pfp];
//
//            if( pfParticlesToShowerKalmanMap.count(pfp) == 0 ){
//                std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\t Warning, no match for a Kalman track for this PFP."<<std::endl;
//                continue;
//            }
//            const art::Ptr<recob::Track> kalman = pfParticlesToShowerKalmanMap[pfp];
//
//            if(kalmanTrackToCaloMap.count(kalman)==0){
//                std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\t Warning, no match for a Calo for this Kalman track."<<std::endl;
//                m_reco_shower_kalman_exists[i_shr]=0;
//                continue;
//            }
//
//            m_reco_shower_kalman_exists[i_shr]=1;
//            const std::vector<art::Ptr<anab::Calorimetry>> calo = kalmanTrackToCaloMap[kalman];
//
//            if(calo.size()!=3){
//                std::cout<<"Singlephoton::AnalyzeKalmanShowerrs\t||\tERROR!! ERROR!!! anab::Calorimetery vector is not of length 3!!! "<<std::endl;
//                continue;
//            }
//            //std::cout<<"index 0 is "<<calo[0]->PlaneID()<<", 1 is "<<calo[1]->PlaneID()<<", 2 is "<<calo[2]->PlaneID()<<std::endl;
//
//            double res_range_lim = m_length_dqdx_box; //4cm 
//
//            for(size_t p=0; p<calo.size();p++){
//    
//                int plane = calo[p]->PlaneID().Plane;
//                if(plane<0 || plane > 3) continue; // Guanqun: plane == 3 is allowed??
//
//                std::vector<double> t_dEdx; //in XX cm only  (4 for now)
//                std::vector<double> t_res;
//            
//
//                for(size_t ix=0; ix<calo[p]->ResidualRange().size(); ++ix){
//
//                    double rr = calo[p]->ResidualRange().back() - calo[p]->ResidualRange()[ix]; 
//                    if(rr <= res_range_lim){
//			// Guanqun: why is here a gains[plane], it's not converting ADC to E?
//                        t_dEdx.push_back(gains[plane]*m_work_function*calo[p]->dQdx()[ix]*1e-6 /m_recombination_factor);
//                        //t_dQdx.push_back(*calo[p]->dQdx()[x]);
//                        t_res.push_back(rr);
//                    }
//                }
//                /*std::cout<<"KAL: plane "<<calo[p]->PlaneID()<<" En: "<<std::endl;
//                for(auto &xx: t_res) std::cout<<xx<<" ";
//                std::cout<<std::endl;
//                //for(auto &xx: t_dEdx) std::cout<<xx<<" ";
//                std::cout<<std::endl<<" dQdx: "<<std::endl;
//                //for(auto &xx: t_dQdx) std::cout<<xx<<" ";
//                std::cout<<std::endl;
//                */
//
//                double tmean = NAN;
//                double tmedian = NAN;
//
//
//                if(t_dEdx.size()>0) tmedian = this->getMedian(t_dEdx);
//                if(t_dEdx.size()>0) tmean = std::accumulate(t_dEdx.begin(), t_dEdx.end(), 0)/((double)t_dEdx.size());
//			std::cout<<"CHECK "<<__LINE__<<" in analyze_Showers.h to see why dEdx does not work"<<std::endl;
//                switch(plane){
//                    case 0:
//                            m_reco_shower_kalman_mean_dEdx_plane0[i_shr] = tmean ;
//                            m_reco_shower_kalman_median_dEdx_plane0[i_shr] = tmedian ;
//                            break;
//                    case 1:
//                            m_reco_shower_kalman_mean_dEdx_plane1[i_shr] = tmean;
//                            m_reco_shower_kalman_median_dEdx_plane1[i_shr] = tmedian;
//                            break;
//                    case 2:
//                            m_reco_shower_kalman_mean_dEdx_plane2[i_shr] = tmean;
//                            m_reco_shower_kalman_median_dEdx_plane2[i_shr] = tmedian;
//                            break;
//                    default:
//                            break;
//                }
//
//
//                const std::vector< anab::Point_t >  kal_pts = calo[p]->XYZ();   
//                double circle = 1.0;//in cm
//                std::vector<double> pts_within;
//                std::vector<double> pts_x;
//
//                for(size_t ix=0; ix<kal_pts.size(); ++ix){
//                    //std::cout<<"KAL: "<<kal_pts[ix].X()<<" "<<kal_pts[ix].Y()<<" "<<kal_pts[ix].Z()<<std::endl;
//                    pts_within.push_back(0);   
//                    pts_x.push_back(calo[p]->ResidualRange().back()-calo[p]->ResidualRange()[ix]);
//
//                    double wire = (double)calcWire(kal_pts[ix].Y(), kal_pts[ix].Z(), plane, m_TPC, m_Cryostat, *geom);
////                    double time = calcTime(kal_pts[ix].X(), plane, m_TPC,m_Cryostat, theDetector);
//					double time = theDetector.ConvertXToTicks(kal_pts[ix].X(), plane, m_TPC,m_Cryostat);
//
//                    //loop over all hits  
//                    for(auto &hit: hitz){
//                        if(plane != hit->View())continue;
//                        double this_w = (double)hit->WireID().Wire;
//                        double this_t = (double)hit->PeakTime();
//                        double dist = sqrt(pow(wire*0.3-this_w*0.3,2)+pow(time/25.0-this_t/25.0,2));
//                        if(dist<=circle) pts_within.back()++;
//                    }
//                    //std::cout<<"KAL "<<ix<<" "<<pts_within.back()<<" "<<calo[p]->ResidualRange().back()-calo[p]->ResidualRange()[ix]<<std::endl;
//                }
//                if(false && pts_x.size()>2){
//                TCanvas *c = new TCanvas();
//                c->cd();
//               
//                TGraph *g  = new TGraph(pts_x.size(), &pts_x[0], &pts_within[0]);
//                g->SetLineColor(kRed);
//                g->SetLineWidth(2);
//                g->Draw("alp");
//                g->SetTitle(("kal_"+std::to_string(plane)+"_"+std::to_string(i_shr)+"_"+std::to_string(m_event_number) +".pdf").c_str());
//                c->SaveAs(("kal_"+std::to_string(plane)+"_"+std::to_string(i_shr)+"_"+std::to_string(m_event_number) +".pdf").c_str(),"pdf");
//                }
//            }
//
//
//            // some kalman averaging
//            double tmp_kal_2 = m_reco_shower_kalman_median_dEdx_plane2[i_shr];
//            double tmp_kal_1 = m_reco_shower_kalman_median_dEdx_plane1[i_shr];
//            double tmp_kal_0 = m_reco_shower_kalman_median_dEdx_plane0[i_shr];
//            double wei_0 = fabs(cos(m_reco_shower_angle_wrt_wires_plane0[i_shr]));
//            double wei_1 = fabs(cos(m_reco_shower_angle_wrt_wires_plane1[i_shr]));
//            double wei_2 = 20.0*fabs(cos(m_reco_shower_angle_wrt_wires_plane2[i_shr]));
//           
//            double thresh = 0.01;
//
//
//            if(tmp_kal_2!=tmp_kal_2 || tmp_kal_2< thresh){
//                tmp_kal_2 = 0;
//                wei_2 = 0.0;
//            }
//            if(tmp_kal_1!=tmp_kal_1 || tmp_kal_1 < thresh){
//                tmp_kal_1 = 0;
//                wei_1 = 0.0;
//            }
//            if(tmp_kal_0!=tmp_kal_0 || tmp_kal_0 < thresh){
//                tmp_kal_0 = 0;
//                wei_0 = 0.0;
//            }
//            double kal_norm = wei_0+wei_1+wei_2;
//
//            if(kal_norm!=0.0){
//                m_reco_shower_kalman_median_dEdx_allplane[i_shr] = (tmp_kal_2*wei_2+tmp_kal_1*wei_1+tmp_kal_0*wei_0)/(kal_norm);
//            }else{
//                m_reco_shower_kalman_median_dEdx_allplane[i_shr] = NAN;
//            }
//
//			std::cout<<"CHECK wei0 "<<wei_0<<" wei1 "<< wei_1<<" wei2 "<<wei_2<<std::endl;
//			std::cout<<i_shr<<" dedx "<<m_reco_shower_kalman_median_dEdx_allplane[i_shr]<<std::endl;
//			sleep(2);
//
//        i_shr++;
//        }
//        return;
    }

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
