namespace single_photon
{

    void SinglePhoton::ClearSecondShowers(){
        m_sss_num_unassociated_hits=0;
        m_sss_num_unassociated_hits_below_threshold=0;
        m_sss_num_associated_hits=0;

        m_sss_num_candidates = 0;

	m_sss_candidate_in_nu_slice.clear();
        m_sss_candidate_num_hits.clear();
        m_sss_candidate_num_wires.clear();
        m_sss_candidate_num_ticks.clear();
        m_sss_candidate_plane.clear();
        m_sss_candidate_PCA.clear();
        m_sss_candidate_mean_ADC.clear();
	m_sss_candidate_ADC_RMS.clear();
        m_sss_candidate_impact_parameter.clear();
        m_sss_candidate_fit_slope.clear();
        m_sss_candidate_veto_score.clear();
        m_sss_candidate_fit_constant.clear();
        m_sss_candidate_mean_tick.clear();
        m_sss_candidate_max_tick.clear();
        m_sss_candidate_min_tick.clear();
        m_sss_candidate_min_wire.clear();
        m_sss_candidate_max_wire.clear();
        m_sss_candidate_mean_wire.clear();
        m_sss_candidate_min_dist.clear();
	m_sss_candidate_wire_tick_based_length.clear();
        m_sss_candidate_energy.clear();
        m_sss_candidate_angle_to_shower.clear();
        m_sss_candidate_closest_neighbour.clear();
        m_sss_candidate_matched.clear();
	m_sss_candidate_matched_energy_fraction_best_plane.clear();
        m_sss_candidate_pdg.clear();
        m_sss_candidate_parent_pdg.clear();
        m_sss_candidate_trackid.clear();
	m_sss_candidate_true_energy.clear();
        m_sss_candidate_overlay_fraction.clear();
        m_sss_candidate_remerge.clear();
    }

    void SinglePhoton::ClearStubs(){
	m_trackstub_num_unassociated_hits = 0; 
        m_trackstub_unassociated_hits_below_threshold = 0; 
        m_trackstub_associated_hits=0; 
        m_trackstub_num_candidates=0; 
	m_trackstub_candidate_in_nu_slice.clear();
        m_trackstub_candidate_num_hits.clear();
        m_trackstub_candidate_num_wires.clear(); 
        m_trackstub_candidate_num_ticks.clear();
        m_trackstub_candidate_plane.clear(); 
        m_trackstub_candidate_PCA.clear();
        m_trackstub_candidate_mean_ADC.clear();
        m_trackstub_candidate_ADC_RMS.clear();
        m_trackstub_candidate_veto_score.clear();
        m_trackstub_candidate_mean_tick.clear();
        m_trackstub_candidate_max_tick.clear();
        m_trackstub_candidate_min_tick.clear();
        m_trackstub_candidate_min_wire.clear();
        m_trackstub_candidate_max_wire.clear();
        m_trackstub_candidate_mean_wire.clear();
        m_trackstub_candidate_min_dist.clear();  
        m_trackstub_candidate_min_impact_parameter_to_shower.clear(); 
        m_trackstub_candidate_min_conversion_dist_to_shower_start.clear();  
        m_trackstub_candidate_min_ioc_to_shower_start.clear();        
        m_trackstub_candidate_ioc_based_length.clear();         
        m_trackstub_candidate_wire_tick_based_length.clear();           
        m_trackstub_candidate_mean_ADC_first_half.clear();              
        m_trackstub_candidate_mean_ADC_second_half.clear();
        m_trackstub_candidate_mean_ADC_first_to_second_ratio.clear(); 
        m_trackstub_candidate_track_angle_wrt_shower_direction.clear();   
        m_trackstub_candidate_linear_fit_chi2.clear();          
        m_trackstub_candidate_energy.clear();
        m_trackstub_candidate_remerge.clear(); 
        m_trackstub_candidate_matched.clear(); 
        m_trackstub_candidate_matched_energy_fraction_best_plane.clear(); 
        m_trackstub_candidate_pdg.clear();   
        m_trackstub_candidate_parent_pdg.clear();
        m_trackstub_candidate_trackid.clear(); 
	m_trackstub_candidate_true_energy.clear();
        m_trackstub_candidate_overlay_fraction.clear(); 

        m_trackstub_num_candidate_groups = 0;                
        m_grouped_trackstub_candidate_indices.clear(); 
        m_trackstub_candidate_group_timeoverlap_fraction.clear();   
    }

    void SinglePhoton::ResizeSecondShowers(size_t size){

    }


    void SinglePhoton::CreateSecondShowerBranches(){
        vertex_tree->Branch("sss_num_unassociated_hits",&m_sss_num_unassociated_hits,"sss_num_unassociated_hits/I");
        vertex_tree->Branch("sss_num_unassociated_hits_below_threshold",&m_sss_num_unassociated_hits_below_threshold,"sss_num_unassociated_hits_below_threshold/I");
        vertex_tree->Branch("sss_num_associated_hits",&m_sss_num_associated_hits,"sss_num_associated_hits/I");

        vertex_tree->Branch("sss_num_candidates",&m_sss_num_candidates,"sss_num_candidates/I");
        vertex_tree->Branch("sss_candidate_veto_score",&m_sss_candidate_veto_score);
	vertex_tree->Branch("sss_candidate_in_nu_slice", &m_sss_candidate_in_nu_slice);
        vertex_tree->Branch("sss_candidate_num_hits",&m_sss_candidate_num_hits);
        vertex_tree->Branch("sss_candidate_num_wires",&m_sss_candidate_num_wires);
        vertex_tree->Branch("sss_candidate_num_ticks",&m_sss_candidate_num_ticks);
        vertex_tree->Branch("sss_candidate_plane",&m_sss_candidate_plane);
        vertex_tree->Branch("sss_candidate_PCA",&m_sss_candidate_PCA);
        vertex_tree->Branch("sss_candidate_mean_ADC",&m_sss_candidate_mean_ADC);
	vertex_tree->Branch("sss_candidate_ADC_RMS", &m_sss_candidate_ADC_RMS);
        vertex_tree->Branch("sss_candidate_impact_parameter",&m_sss_candidate_impact_parameter); 
        vertex_tree->Branch("sss_candidate_fit_slope",&m_sss_candidate_fit_slope);
        vertex_tree->Branch("sss_candidate_fit_constant",&m_sss_candidate_fit_constant);
        vertex_tree->Branch("sss_candidate_mean_tick",&m_sss_candidate_mean_tick);
        vertex_tree->Branch("sss_candidate_max_tick",&m_sss_candidate_max_tick);
        vertex_tree->Branch("sss_candidate_min_tick",&m_sss_candidate_min_tick);
        vertex_tree->Branch("sss_candidate_mean_wire",&m_sss_candidate_mean_wire);
        vertex_tree->Branch("sss_candidate_max_wire",&m_sss_candidate_max_wire);
        vertex_tree->Branch("sss_candidate_min_wire",&m_sss_candidate_min_wire);
        vertex_tree->Branch("sss_candidate_min_dist",&m_sss_candidate_min_dist);
	vertex_tree->Branch("sss_candidate_wire_tick_based_length", &m_sss_candidate_wire_tick_based_length);
        vertex_tree->Branch("sss_candidate_energy",&m_sss_candidate_energy);
        vertex_tree->Branch("sss_candidate_angle_to_shower",&m_sss_candidate_angle_to_shower);
        vertex_tree->Branch("sss_candidate_closest_neighbour",&m_sss_candidate_closest_neighbour);
        vertex_tree->Branch("sss_candidate_remerge",&m_sss_candidate_remerge);

        vertex_tree->Branch("sss_candidate_matched",&m_sss_candidate_matched);
        vertex_tree->Branch("sss_candidate_pdg",&m_sss_candidate_pdg);
        vertex_tree->Branch("sss_candidate_parent_pdg",&m_sss_candidate_parent_pdg);
        vertex_tree->Branch("sss_candidate_trackid",&m_sss_candidate_trackid);
	vertex_tree->Branch("sss_candidate_true_energy", &m_sss_candidate_true_energy);
        vertex_tree->Branch("sss_candidate_overlay_fraction",&m_sss_candidate_overlay_fraction);
	vertex_tree->Branch("sss_candidate_matched_energy_fraction_best_plane", &m_sss_candidate_matched_energy_fraction_best_plane);


        vertex_tree->Branch("sss3d_ioc_ranked_en",&m_sss3d_ioc_ranked_en);
        vertex_tree->Branch("sss3d_ioc_ranked_conv",&m_sss3d_ioc_ranked_conv);
        vertex_tree->Branch("sss3d_ioc_ranked_invar",&m_sss3d_ioc_ranked_invar);
        vertex_tree->Branch("sss3d_ioc_ranked_implied_invar",&m_sss3d_ioc_ranked_implied_invar);
        vertex_tree->Branch("sss3d_ioc_ranked_ioc",&m_sss3d_ioc_ranked_ioc);
        vertex_tree->Branch("sss3d_ioc_ranked_opang",&m_sss3d_ioc_ranked_opang);
        vertex_tree->Branch("sss3d_ioc_ranked_implied_opang",&m_sss3d_ioc_ranked_implied_opang);
        vertex_tree->Branch("sss3d_ioc_ranked_id",&m_sss3d_ioc_ranked_id);

        vertex_tree->Branch("sss3d_invar_ranked_en",&m_sss3d_invar_ranked_en);
        vertex_tree->Branch("sss3d_invar_ranked_conv",&m_sss3d_invar_ranked_conv);
        vertex_tree->Branch("sss3d_invar_ranked_invar",&m_sss3d_invar_ranked_invar);
        vertex_tree->Branch("sss3d_invar_ranked_implied_invar",&m_sss3d_invar_ranked_implied_invar);
        vertex_tree->Branch("sss3d_invar_ranked_ioc",&m_sss3d_invar_ranked_ioc);
        vertex_tree->Branch("sss3d_invar_ranked_opang",&m_sss3d_invar_ranked_opang);
        vertex_tree->Branch("sss3d_invar_ranked_implied_opang",&m_sss3d_invar_ranked_implied_opang);
        vertex_tree->Branch("sss3d_invar_ranked_id",&m_sss3d_invar_ranked_id);


        vertex_tree->Branch("sss2d_ioc_ranked_en",&m_sss2d_ioc_ranked_en);
        vertex_tree->Branch("sss2d_ioc_ranked_conv",&m_sss2d_ioc_ranked_conv);
        vertex_tree->Branch("sss2d_ioc_ranked_ioc",&m_sss2d_ioc_ranked_ioc);
        vertex_tree->Branch("sss2d_ioc_ranked_pca",&m_sss2d_ioc_ranked_pca);
        vertex_tree->Branch("sss2d_ioc_ranked_invar",&m_sss2d_ioc_ranked_invar);
        vertex_tree->Branch("sss2d_ioc_ranked_angle_to_shower",&m_sss2d_ioc_ranked_angle_to_shower);
        vertex_tree->Branch("sss2d_ioc_ranked_num_planes",&m_sss2d_ioc_ranked_num_planes);

        vertex_tree->Branch("sss2d_invar_ranked_en",&m_sss2d_invar_ranked_en);
        vertex_tree->Branch("sss2d_invar_ranked_conv",&m_sss2d_invar_ranked_conv);
        vertex_tree->Branch("sss2d_invar_ranked_ioc",&m_sss2d_invar_ranked_ioc);
        vertex_tree->Branch("sss2d_invar_ranked_pca",&m_sss2d_invar_ranked_pca);
        vertex_tree->Branch("sss2d_invar_ranked_invar",&m_sss2d_invar_ranked_invar);
        vertex_tree->Branch("sss2d_invar_ranked_angle_to_shower",&m_sss2d_invar_ranked_angle_to_shower);
        vertex_tree->Branch("sss2d_invar_ranked_num_planes",&m_sss2d_invar_ranked_num_planes);

        vertex_tree->Branch("sss2d_conv_ranked_en",&m_sss2d_conv_ranked_en);
        vertex_tree->Branch("sss2d_conv_ranked_conv",&m_sss2d_conv_ranked_conv);
        vertex_tree->Branch("sss2d_conv_ranked_ioc",&m_sss2d_conv_ranked_ioc);
        vertex_tree->Branch("sss2d_conv_ranked_pca",&m_sss2d_conv_ranked_pca);
        vertex_tree->Branch("sss2d_conv_ranked_invar",&m_sss2d_conv_ranked_invar);
        vertex_tree->Branch("sss2d_conv_ranked_angle_to_shower",&m_sss2d_conv_ranked_angle_to_shower);
        vertex_tree->Branch("sss2d_conv_ranked_num_planes",&m_sss2d_conv_ranked_num_planes);

    }

    void SinglePhoton::CreateStubBranches(){

        vertex_tree->Branch("trackstub_num_unassociated_hits",&m_trackstub_num_unassociated_hits,"trackstub_num_unassociated_hits/I");
        vertex_tree->Branch("trackstub_unassociated_hits_below_threshold",&m_trackstub_unassociated_hits_below_threshold,"trackstub_unassociated_hits_below_threshold/I");
        vertex_tree->Branch("trackstub_associated_hits",&m_trackstub_associated_hits,"trackstub_associated_hits/I");
	vertex_tree->Branch("trackstub_num_candidates", &m_trackstub_num_candidates, "trackstub_num_candidates/I");
	vertex_tree->Branch("trackstub_candidate_in_nu_slice", &m_trackstub_candidate_in_nu_slice);
	vertex_tree->Branch("trackstub_candidate_num_hits", &m_trackstub_candidate_num_hits);
	vertex_tree->Branch("trackstub_candidate_num_wires", &m_trackstub_candidate_num_wires);
	vertex_tree->Branch("trackstub_candidate_num_ticks", &m_trackstub_candidate_num_ticks);
	vertex_tree->Branch("trackstub_candidate_plane", &m_trackstub_candidate_plane);
	vertex_tree->Branch("trackstub_candidate_PCA", &m_trackstub_candidate_PCA);
	vertex_tree->Branch("trackstub_candidate_mean_ADC", &m_trackstub_candidate_mean_ADC);
	vertex_tree->Branch("trackstub_candidate_ADC_RMS", &m_trackstub_candidate_ADC_RMS);
	vertex_tree->Branch("trackstub_candidate_veto_score", &m_trackstub_candidate_veto_score);
	vertex_tree->Branch("trackstub_candidate_mean_tick", &m_trackstub_candidate_mean_tick);
	vertex_tree->Branch("trackstub_candidate_max_tick", &m_trackstub_candidate_max_tick);
	vertex_tree->Branch("trackstub_candidate_min_tick", &m_trackstub_candidate_min_tick);
	vertex_tree->Branch("trackstub_candidate_min_wire", &m_trackstub_candidate_min_wire);
	vertex_tree->Branch("trackstub_candidate_max_wire", &m_trackstub_candidate_max_wire);
	vertex_tree->Branch("trackstub_candidate_mean_wire", &m_trackstub_candidate_mean_wire);
	vertex_tree->Branch("trackstub_candidate_min_dist", &m_trackstub_candidate_min_dist);
	vertex_tree->Branch("trackstub_candidate_min_impact_parameter_to_shower", &m_trackstub_candidate_min_impact_parameter_to_shower);
	vertex_tree->Branch("trackstub_candidate_min_conversion_dist_to_shower_start", &m_trackstub_candidate_min_conversion_dist_to_shower_start);
	vertex_tree->Branch("trackstub_candidate_min_ioc_to_shower_start", &m_trackstub_candidate_min_ioc_to_shower_start);
	vertex_tree->Branch("trackstub_candidate_ioc_based_length", &m_trackstub_candidate_ioc_based_length);
	vertex_tree->Branch("trackstub_candidate_wire_tick_based_length", &m_trackstub_candidate_wire_tick_based_length);
	vertex_tree->Branch("trackstub_candidate_mean_ADC_first_half", &m_trackstub_candidate_mean_ADC_first_half);
	vertex_tree->Branch("trackstub_candidate_mean_ADC_second_half", &m_trackstub_candidate_mean_ADC_second_half);
	vertex_tree->Branch("trackstub_candidate_mean_ADC_first_to_second_ratio", &m_trackstub_candidate_mean_ADC_first_to_second_ratio);
	vertex_tree->Branch("trackstub_candidate_track_angle_wrt_shower_direction", &m_trackstub_candidate_track_angle_wrt_shower_direction);
	vertex_tree->Branch("trackstub_candidate_linear_fit_chi2", &m_trackstub_candidate_linear_fit_chi2);
	vertex_tree->Branch("trackstub_candidate_energy", &m_trackstub_candidate_energy);
	vertex_tree->Branch("trackstub_candidate_remerge", &m_trackstub_candidate_remerge);
	vertex_tree->Branch("trackstub_candidate_matched", &m_trackstub_candidate_matched);
	vertex_tree->Branch("trackstub_candidate_matched_energy_fraction_best_plane", &m_trackstub_candidate_matched_energy_fraction_best_plane);
	vertex_tree->Branch("trackstub_candidate_pdg", &m_trackstub_candidate_pdg);
	vertex_tree->Branch("trackstub_candidate_parent_pdg", &m_trackstub_candidate_parent_pdg);
	vertex_tree->Branch("trackstub_candidate_trackid", &m_trackstub_candidate_trackid);
	vertex_tree->Branch("trackstub_candidate_true_energy", &m_trackstub_candidate_true_energy);
	vertex_tree->Branch("trackstub_candidate_overlay_fraction", &m_trackstub_candidate_overlay_fraction);


	vertex_tree->Branch("trackstub_num_candidate_groups", &m_trackstub_num_candidate_groups, "trackstub_num_candidate_groups/I");
	vertex_tree->Branch("grouped_trackstub_candidate_indices", &m_grouped_trackstub_candidate_indices);
	vertex_tree->Branch("trackstub_candidate_group_timeoverlap_fraction", &m_trackstub_candidate_group_timeoverlap_fraction);
    
    }





    TGraph* SinglePhoton::GetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts){

        std::vector<double>t_wire;
        std::vector<double>t_tick;
        // std::vector<double>t_dist;

        std::vector<double>all_wire; // wire of all hits
        std::vector<double>all_tick;
        std::vector<double>all_dist; // distance to vertex of all hits


        for(size_t h = 0; h< hitz.size(); h++){
            auto hit = hitz[h];
            double h_wire = (double)hit->WireID().Wire;
            double h_tick = (double)hit->PeakTime();

            double dd =sqrt(pow(h_wire*0.3-vertex_wire*0.3,2)+pow(h_tick/25.0- vertex_tick/25.0,2));
            all_wire.push_back(h_wire);   
            all_tick.push_back(h_tick);   
            all_dist.push_back(dd);
        }

        std::vector<size_t> sorted_in = sort_indexes(all_dist); // index of all dist in descending order
        size_t max_e = std::min((size_t)Npts,hitz.size());

        for(size_t i =0; i<max_e; i++){
            t_wire.push_back(all_wire[sorted_in[hitz.size()-1-i]]);
            t_tick.push_back(all_tick[sorted_in[hitz.size()-1-i]]);
        }

        return new TGraph(t_wire.size(),&t_wire[0],&t_tick[0]);
    }

    sss_score SinglePhoton::ScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower){
        sss_score score(p,cl);
        score.n_hits = hits.size();

        std::vector<double> t_wires;
        std::vector<double> t_ticks;

        // 
        int n_min_ticks = 4;
        int n_min_wires = 3;
        double n_max_pca = 0.9999;

        score.pass = true;

        // ************* Some simple metrics relative to study point (usually vertex) ***************
        // this can be moved to inclass initializer
        score.max_dist_tick = 0;
        score.min_dist_tick = 1e10;
        score.mean_dist_tick = 0;

        score.max_dist_wire = 0;
        score.min_dist_wire = 1e10;
        score.mean_dist_wire = 0;

        score.max_dist = 0;
        score.min_dist = 1e10;
        score.mean_dist = 0;

        score.mean_tick =0;
        score.max_tick =0;
        score.min_tick =1e10;

        score.mean_wire =0;
        score.max_wire =0;
        score.min_wire =1e10;

        score.n_wires = 0;
        score.n_ticks = 0;

        score.impact_parameter = -99;

        score.close_tick = -99;
        score.close_wire = -99;

        std::map<int,bool> wire_count;
        std::map<int,bool> tick_count;

        for(auto &h: hits){
            double h_tick = (double)h->PeakTime();
            double h_wire = (double)h->WireID().Wire;

            score.mean_wire += h_wire;
            score.mean_tick += h_tick;

            score.max_wire = std::max(score.max_wire, h_wire);
            score.min_wire = std::min(score.min_wire, h_wire);

            score.max_tick = std::max(score.max_tick, h_tick);
            score.min_tick = std::min(score.min_tick, h_tick);

            score.max_dist_tick = std::max(score.max_dist_tick, fabs(h_tick-vertex_tick));
            score.min_dist_tick = std::min(score.min_dist_tick, fabs(h_tick-vertex_tick));

            score.max_dist_wire = std::max(score.max_dist_wire, fabs(h_wire-vertex_wire));
            score.min_dist_wire = std::min(score.min_dist_wire, fabs(h_wire-vertex_wire));

            score.mean_dist_tick += fabs(h_tick-vertex_tick);
            score.mean_dist_wire += fabs(h_wire-vertex_wire);

            //wierd dits
            //requires that hit in hits has to be on the same plane as vertex_wire.
            double dd =sqrt(pow(h_wire*0.3-vertex_wire*0.3,2)+pow(h_tick/25.0- vertex_tick/25.0,2));
            score.mean_dist += dd;
            if(dd< score.min_dist){
                score.close_wire = h_wire;
                score.close_tick = h_tick;
            }

            score.max_dist = std::max(dd,score.max_dist);
            score.min_dist = std::min(dd,score.min_dist);


            t_wires.push_back(h_wire);
            t_ticks.push_back(h_tick);

            if(wire_count.count((int)h_wire)<1){
                wire_count[((int)h_wire)] = true;
                score.n_wires++;
            }
            if(tick_count.count((int)h_tick)<1){
                tick_count[((int)h_tick)] = true;
                score.n_ticks++;
            }

        }

        //            TGraph * g_pts = new TGraph(t_wires.size(),&t_ticks[0],&t_wires[0]);

        score.mean_tick = score.mean_tick/(double)score.n_hits;
        score.mean_wire = score.mean_wire/(double)score.n_hits;

        score.mean_dist = score.mean_dist/(double)score.n_hits;

        score.mean_dist_tick = score.mean_dist_tick/(double)score.n_hits;
        score.mean_dist_wire = score.mean_dist_wire/(double)score.n_hits;

        // **************** Metrics of Pointing: Does this cluster "point" back to the vertex? *************************
        // **************** First off, PCA

        TPrincipal* principal = new TPrincipal(2,"D");
        double mod_wire = 1.0;
        double mod_tick = 1.0;

        for(int i = 0; i < score.n_hits; i++){
            std::vector<double> tmp_pts = {t_wires[i]*mod_wire, t_ticks[i]/mod_tick};
            principal->AddRow(&tmp_pts[0]);
        }
        principal->MakePrincipals();
        //principal->Print();

        TVectorD * eigenval = (TVectorD*) principal->GetEigenValues();
        //TMatrixD * eigenvec = (TMatrixD*) principal->GetEigenVectors();
        TMatrixD * covar = (TMatrixD*) principal->GetCovarianceMatrix();

        score.pca_0 = (*eigenval)(0);
        score.pca_1 = (*eigenval)(1);

        //(*eigenvec).Print();
        //(*covar).Print();
        //std::cout<<"SinglePhoton::SSS\t||\tEigen: "<<score.pca_0<<" "<<score.pca_1<<std::endl;

        score.pca_theta = atan((*covar)[0][0]/(*covar)[0][1])*180.0/3.14159;


        double slope = ((*covar)[0][0]/(*covar)[0][1]);
        double c = score.mean_tick*mod_wire - slope*score.mean_wire/mod_tick;
        score.impact_parameter = fabs(slope*vertex_wire*mod_wire +vertex_tick/mod_tick+c)/sqrt(slope*slope+1.0*1.0);


        if(score.n_wires < n_min_wires || score.n_ticks < n_min_ticks || score.pca_0 >= n_max_pca){
            score.pass = false;
        }




        delete principal;

        return score;
    }

    int SinglePhoton::CompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,
            const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      const   std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,                    double eps){


        for(size_t s =0; s< showers.size(); s++){
            art::Ptr<recob::Shower> shower = showers[s];
            art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(shower);
            std::vector<art::Ptr<recob::Hit>> showerhits = pfParticleToHitsMap.at(pfp);

            bool in_primary_shower = false;
            for(size_t h = 0; h< hitz.size(); h++){
                auto hit = hitz[h];
                double h_wire = (double)hit->WireID().Wire;
                double h_tick = (double)hit->PeakTime();


                for(auto &sh: showerhits){

                    if(sh->View() != hit->View()) continue;

                    double sh_wire = (double)sh->WireID().Wire;
                    double sh_tick = (double)sh->PeakTime();


                    double dist = sqrt(pow(sh_wire*0.3-h_wire*0.3,2)+pow(sh_tick/25.0-h_tick/25.0,2));

                    if(dist<=eps){
                        in_primary_shower = true;
                        return (int)s;
                    }

                }

            }

            if(in_primary_shower){
                return (int)s;
            }
        }


        return -1;
    }



    std::vector<double>SinglePhoton::SecondShowerMatching(
	    std::vector<art::Ptr<recob::Hit>>& hitz,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
            std::map< int ,art::Ptr<simb::MCParticle>>  & MCParticleToTrackIdMap
	){


        std::vector<double> ans; //matched,pdg,parentpdg,trkid


        std::vector<double> vec_fraction_matched;
        std::map<std::string,bool> map_is_shower_process = {{"compt",true},{"FastScintillation",true},{"eBrem",true},{"phot",true},{"eIoni",true},{"conv",true},{"annihil",true}};
        bool reco_verbose = false;

        std::unordered_map<int,double> objide; //map between the MCParticle track ID and the backtracker energy

        //energy for an MCParticle that comprises the most energy when sum over associated hits in PFP
        //total energy of the reco PFP taken from the sum of the hits associated to an MCParticle
        double maxe=-1, tote=0;                

        std::vector<double> total_energy_on_plane = {0.0,0.0,0.0};
        art::Ptr<simb::MCParticle> best_matched_mcparticle; //pointer for the MCParticle match we will calculate
        std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the cluster
        std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

        int n_associated_mcparticle_hits = 0;
        int n_not_associated_hits = 0;

        //this is the vector that will store the associated MC paritcles, as well as a MAP to the amount of energy associated
        std::vector<art::Ptr<simb::MCParticle>> asso_mcparticles_vec;
        std::map<art::Ptr<simb::MCParticle>, std::vector<double>> map_asso_mcparticles_energy;
        bool found_a_match = false;

        //loop only over hits associated to this reco PFP
        for(size_t i_h=0; i_h < hitz.size(); ++i_h){
            int which_plane = (int)hitz[i_h]->View();
            particle_vec.clear(); match_vec.clear(); //only store per hit
            //for the hit, fill the backtracker info 
            mcparticles_per_hit.get(hitz[i_h].key(), particle_vec, match_vec);

            //if there is an MCParticle associated to this hit
            if(particle_vec.size()>0) n_associated_mcparticle_hits++;
            if(particle_vec.size()==0) n_not_associated_hits++;

            //for each MCParticle associated with this hit
            for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
                //add the energy of the back tracked hit for this MCParticle to the track id for the MCParticle in the map
                objide[ particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy; //store energy per track id

                //if the id isn't already in the map, store it in the vector of all associated MCParticles
                if(std::find(asso_mcparticles_vec.begin(), asso_mcparticles_vec.end(),  particle_vec[i_p]) == asso_mcparticles_vec.end()){
                    asso_mcparticles_vec.push_back(particle_vec[i_p]);
                    map_asso_mcparticles_energy[particle_vec[i_p]] = {0.0,0.0,0.0};
                    map_asso_mcparticles_energy[particle_vec[i_p]][which_plane] =  match_vec[i_p]->energy;
                }else{
                    map_asso_mcparticles_energy[particle_vec[i_p]][which_plane] += match_vec[i_p]->energy;
                }
                //add the energy of the back tracked hit to the total energy for the PFP
                tote += match_vec[i_p]->energy; //calculate total energy deposited
                total_energy_on_plane[which_plane]+=match_vec[i_p]->energy;

                //want the MCParticle with the max total energy summed from the back tracker hit energy from hits in PFP
                //TODO: this part will change once the parts below are fully implemented
                if( objide[ particle_vec[i_p]->TrackId()] > maxe ){ //keep track of maximum
                    maxe = objide[ particle_vec[i_p]->TrackId() ];
                    best_matched_mcparticle = particle_vec[i_p]; //we will now define the best match as a source MCP rather than the max single energy contributor 
                    found_a_match = true;//will be false for showers from overlay
                }
            }//end loop over particles per hit

        } // end loop over hit
        if(found_a_match){
            std::cout<<"Found a match!"<<std::endl;
        }
        double fraction_num_hits_overlay = (double)n_not_associated_hits/(double)hitz.size();

        //            if(reco_verbose)std::cout << "SinglePhoton::recoMC()\t||\t On Object "<<i<<". The number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       
        //          if(reco_verbose) std::cout<<"SinglePhoton::recoMC()\t||\t the fraction of hits from overlay is is "<<fraction_num_hits_overlay<<" ("<<n_not_associated_hits<<"/"<<obj_hits_ptrs.size()<<")"<<std::endl;


        if(n_associated_mcparticle_hits == 0){
            //This will only occur if the whole recob::PFParticle is PURELY associated with an overlay object
            found_a_match =false;
            //Here we will fill every sim_shower_XXX variable with -999 or something like that 
            return {0,0,0,0,0,0,0};
        }//

        /*
         *
         * Loop over each MCParticle associated to the reco shower to find the source particle
         *
         */

        std::map<int, art::Ptr<simb::MCParticle>> mother_MCP_map; //map between MCP track id and the source MCP

        std::vector<art::Ptr<simb::MCParticle>> marks_mother_vector;
        std::map<art::Ptr<simb::MCParticle>, std::vector<double>> marks_mother_energy_fraction_map;

        int this_mcp_id = -1; //the track id for the current MCP in parent tree
        int last_mcp_id = -1; //the track id for the previous MCP in parent tree
        int i_mcp = 0;

        int num_bt_mothers =0;

        //reco_verbose = false;
        //for each MCP that's associated to the reco shower
        for(auto mcp:asso_mcparticles_vec){

            if(reco_verbose) std::cout<<"-----------------------------Start L1 Loop --------------------------------------------------"<<std::endl;
            //                if(reco_verbose) std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<") Start by Looking at an MCP with pdg code "<<mcp->PdgCode()<<" and status code "<<mcp->StatusCode()<<" TrackID: "<<mcp->TrackId()<<std::endl;
            //              if(reco_verbose) std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<") This MCP gave "<<   map_asso_mcparticles_energy[mcp][0] <<" | "<<map_asso_mcparticles_energy[mcp][1]<<" | "<<map_asso_mcparticles_energy[mcp][2]<<" energy to the recob::Object on each plane"<<std::endl;
            //                std::cout<<"L1: the mother of this MCP is track id "<<mcp->Mother()<<" and there are "<<mcp->NumberDaughters()<<" daughters"<<std::endl;

            //get the track ID for the current MCP
            this_mcp_id = mcp->TrackId();
            last_mcp_id = this_mcp_id;//initialize the previous one

            //while the track id is valid, move up the parent tree for the MCP that contributes to the reco shower
            //currently it keeps going until it hits the top of the interaction chain, but this is likely too far
            //to do a proper match you need to check for different cases and stop once one is fulfilled
            while(this_mcp_id >= 0 ){                  
                art::Ptr<simb::MCParticle> this_mcp = MCParticleToTrackIdMap[this_mcp_id];//get the MCP associated to the track ID
                // std::cout<<"going up the tree got mother particle"<<std::endl;

                //check if it's a valid MCP
                if (this_mcp.isNull()){
                    //                       if(reco_verbose)   std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<")  null pointer at id "<<this_mcp_id<<std::endl;
                    this_mcp_id = last_mcp_id; //if invalid, move back a level to the previous MCP in parent tree and break the loop
                    break;
                }

                //If primary particle will have process "primary"
                //                    if(reco_verbose)    std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<")  going up the tree at an MCP with track id  "<<this_mcp_id<<", pdg code "<<this_mcp->PdgCode()<<", and status code "<<this_mcp->StatusCode()<<" and Mother: "<<this_mcp->Mother()<<" Process: "<<this_mcp->Process()<<" EndProcess: "<<this_mcp->EndProcess()<<std::endl;

                //if it is a valid particle, iterate forward to the mother
                last_mcp_id = this_mcp_id;
                this_mcp_id =  this_mcp->Mother();

                //Check to see if this MCP was created in a "showery" process
                if(map_is_shower_process.count(this_mcp->Process()) > 0){
                    //if it was, keep going, 

                }else if(this_mcp->Process()=="primary"){
                    //if its primary, great! Note it.
                    if(reco_verbose)  std::cout<<"L1: Backtracked to primary! breaking"<<std::endl;
                    this_mcp_id = last_mcp_id; //if invalid, move back a level to the previous MCP in parent tree and break the loop
                    break;
                }else{
                    if(reco_verbose) std::cout<<"L1: Backtracked to a particle created in "<<this_mcp->EndProcess()<<"! breaking"<<std::endl;
                    this_mcp_id = last_mcp_id; //if invalid, move back a level to the previous MCP in parent tree and break the loop
                    break;
                }
            }

            //if the MCP at the top of the interaction chain has a valid track id store this in the mother map
            if (this_mcp_id >= 0){

                mother_MCP_map[this_mcp_id] = MCParticleToTrackIdMap[this_mcp_id];//putting it in a map allows for multiple contributing MCP's in the reco shower to have the same mother MCP

                bool is_old = false;

                for(size_t k=0; k< marks_mother_vector.size(); k++){
                    //if its in it before, just run with it
                    if(marks_mother_vector[k]==MCParticleToTrackIdMap[this_mcp_id]){
                        marks_mother_energy_fraction_map[marks_mother_vector[k]][0] += map_asso_mcparticles_energy[mcp][0];
                        marks_mother_energy_fraction_map[marks_mother_vector[k]][1] += map_asso_mcparticles_energy[mcp][1];
                        marks_mother_energy_fraction_map[marks_mother_vector[k]][2] += map_asso_mcparticles_energy[mcp][2];
                        is_old = true;
                        break;
                    }
                }
                if(is_old==false){
                    marks_mother_vector.push_back(MCParticleToTrackIdMap[this_mcp_id]);
                    marks_mother_energy_fraction_map[marks_mother_vector.back()] = {0.0,0.0,0.0};
                    marks_mother_energy_fraction_map[marks_mother_vector.back()][0] =  map_asso_mcparticles_energy[mcp][0];
                    marks_mother_energy_fraction_map[marks_mother_vector.back()][1] =  map_asso_mcparticles_energy[mcp][1];
                    marks_mother_energy_fraction_map[marks_mother_vector.back()][2] =  map_asso_mcparticles_energy[mcp][2];
                }


                num_bt_mothers++;
            } else{
                if(reco_verbose)  std::cout<<"L1: error, the mother mother id was "<<this_mcp_id <<std::endl;

            }

            if(reco_verbose)  std::cout<<"-----------------------------End L1 Loop --------------------------------------------------"<<std::endl;
            i_mcp++;
        }//for each MCParticle that's associated to a the recob::Shower

        //reco_verbose = true;
        //there should be at least 1 mother MCP
        if(reco_verbose)           std::cout<<"SinglePhoton::recoMC()\t||\t the number of source mother particles is "<<mother_MCP_map.size()<<" of which : "<<marks_mother_vector.size()<<" are unique!"<<std::endl;

        if(reco_verbose)       std::cout<<"---------------------------- L2-------------------------------"<<std::endl;

        double best_mother_index = 0;
        double best_mother_energy = -9999;
        int best_mother_plane = -99;

        for(size_t p=0; p< marks_mother_vector.size(); p++){
            art::Ptr<simb::MCParticle> mother = marks_mother_vector[p];
            std::vector<double> mother_energy_recod = marks_mother_energy_fraction_map[mother];
            if(reco_verbose)    std::cout<<"L2: Mother candidate "<<p<<" TrackID "<<mother->TrackId()<<" Process: "<<mother->Process()<<" EndProcess: "<<mother->EndProcess()<<std::endl;
            if(reco_verbose)   std::cout<<"L2: Mother candidate "<<p<<" Energy "<<mother->E()<<" Reco'd Energy: "<<mother_energy_recod[0]<<" | "<<mother_energy_recod[1]<<" | "<<mother_energy_recod[2]<<" Fraction: ("<<mother_energy_recod[0]/(1000*mother->E())*100.0<<"% , "<<mother_energy_recod[1]/(1000*mother->E())*100.0<<"% , "<<mother_energy_recod[2]/(1000*mother->E())*100.0<<"% )"<<std::endl;

            if( mother_energy_recod[0] > best_mother_energy){
                best_mother_index = p;
                best_mother_energy = mother_energy_recod[0];
                best_mother_plane = 0;
            }

            if( mother_energy_recod[1] > best_mother_energy){
                best_mother_index = p;
                best_mother_energy = mother_energy_recod[1];
                best_mother_plane = 1;
            }

            if( mother_energy_recod[2] > best_mother_energy){
                best_mother_index = p;
                best_mother_energy = mother_energy_recod[2];
                best_mother_plane = 2;
            }

        }

        if(marks_mother_vector.size()!=0){
            //if(reco_verbose)  std::cout<<"SinglePhoton::recoMC()\t||\t The `BEST` mother is a "<<marks_mother_vector[best_mother_index]->PdgCode()<<" at "<<best_mother_index<<" on plane: "<<best_mother_plane<<std::endl;
            std::cout<<"SinglePhoton::recoMC()\t||\t The `BEST` mother is a "<<marks_mother_vector[best_mother_index]->PdgCode()<<" at "<<best_mother_index<<" on plane: "<<best_mother_plane<<std::endl;
            for(int l=0; l<3; l++){
                std::cout<<"SinglePhoton::recoMC()\t||\t It represents "<<marks_mother_energy_fraction_map[marks_mother_vector[best_mother_index]][l]/total_energy_on_plane[l]*100.0<<"% of the energy on plane: "<<l<<" which is "<<total_energy_on_plane[l] <<std::endl;
            }
        }


        if(reco_verbose) std::cout<<"---------------------------- L2-------------------------------"<<std::endl;
        const art::Ptr<simb::MCParticle> match = marks_mother_vector[best_mother_index];

        std::vector<double> corrected_vertex(3), corrected_start(3);
        this->spacecharge_correction(match, corrected_vertex);


        art::Ptr<simb::MCParticle> match_mother = MCParticleToTrackIdMap[match->Mother()];
        int par_pdg = -1;
        if (match_mother.isNull()){
            par_pdg = -1;

        }else{
            par_pdg = match_mother->PdgCode();
        }

        ans = {1,(double)match->PdgCode(), (double)par_pdg, (double)match->TrackId(), match->E(), fraction_num_hits_overlay, best_mother_energy/total_energy_on_plane.at(best_mother_plane)};

        return ans;
    }//end sss matching;






    //************************************************ Shower Search Slice Second SSS3D ********** /

    void SinglePhoton::ClearSecondShowers3D(){

        m_sss3d_num_showers = 0;
        m_sss3d_shower_start_x.clear();
        m_sss3d_shower_start_y.clear();
        m_sss3d_shower_start_z.clear();
        m_sss3d_shower_dir_x.clear();
        m_sss3d_shower_dir_y.clear();
        m_sss3d_shower_dir_z.clear();
        m_sss3d_shower_length.clear();
        m_sss3d_shower_conversion_dist.clear();
        m_sss3d_shower_invariant_mass.clear();
        m_sss3d_shower_implied_invariant_mass.clear();
        m_sss3d_shower_impact_parameter.clear();
        m_sss3d_shower_energy_max.clear();
        m_sss3d_shower_score.clear();
        m_sss3d_slice_nu.clear();
        m_sss3d_slice_clear_cosmic.clear();
        m_sss3d_shower_ioc_ratio.clear();
    }


    void SinglePhoton::CreateSecondShowerBranches3D(){
        vertex_tree->Branch("sss3d_num_showers",&m_sss3d_num_showers,"sss3d_num_showers/I");

        vertex_tree->Branch("sss3d_shower_start_x",&m_sss3d_shower_start_x);
        vertex_tree->Branch("sss3d_shower_start_y",&m_sss3d_shower_start_y);
        vertex_tree->Branch("sss3d_shower_start_z",&m_sss3d_shower_start_z);
        vertex_tree->Branch("sss3d_shower_dir_x",&m_sss3d_shower_dir_x);
        vertex_tree->Branch("sss3d_shower_dir_y",&m_sss3d_shower_dir_y);
        vertex_tree->Branch("sss3d_shower_dir_z",&m_sss3d_shower_dir_z);

        vertex_tree->Branch("sss3d_shower_length",&m_sss3d_shower_length);
        vertex_tree->Branch("sss3d_shower_conversion_dist",&m_sss3d_shower_conversion_dist);
        vertex_tree->Branch("sss3d_shower_invariant_mass",&m_sss3d_shower_invariant_mass);
        vertex_tree->Branch("sss3d_shower_implied_invariant_mass",&m_sss3d_shower_implied_invariant_mass);
        vertex_tree->Branch("sss3d_shower_impact_parameter",&m_sss3d_shower_impact_parameter);
        vertex_tree->Branch("sss3d_shower_ioc_ratio",&m_sss3d_shower_ioc_ratio);
        vertex_tree->Branch("sss3d_shower_energy_max",&m_sss3d_shower_energy_max);
        vertex_tree->Branch("sss3d_shower_score",&m_sss3d_shower_score);
        //vertex_tree->Branch("sss3d_slice_nu",&m_sss3d_slice_nu);
        //vertex_tree->Branch("sss3d_slice_clear_cosmic",&m_sss3d_slice_clear_cosmic);
    }


    void SinglePhoton::SecondShowerSearch3D(std::vector<art::Ptr<recob::Shower>> & showers,std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & NormalShowerToPFParticleMap,  std::vector<art::Ptr<recob::Track>> & tracks, std::map<art::Ptr<recob::Track>,  art::Ptr<recob::PFParticle>> & NormalTrackToPFParticleMap, art::Event const & evt ){

        std::string sss3dlabel = "pandoraShower";//"pandoraAllOutcomesShower"
//Keng        std::string sss3dlabel = "allShr";//"pandoraAllOutcomesShower"
        double max_conv_dist = 80.0;

        art::ValidHandle<std::vector<recob::Shower>> const & allShowerHandle  = evt.getValidHandle<std::vector<recob::Shower>>(sss3dlabel);
        std::vector<art::Ptr<recob::Shower>> allShowerVector;
        art::fill_ptr_vector(allShowerVector,allShowerHandle);
        std::cout<<"We have "<<showers.size()<<" showers in primary slice and "<<allShowerVector.size()<<" in full event."<<std::endl;

        art::FindManyP<recob::Hit> hits_per_shower(allShowerHandle, evt, sss3dlabel);
        std::map<art::Ptr<recob::Shower>, std::vector<art::Ptr<recob::Hit>> > showerToHitsMap;
        for(size_t i=0; i< allShowerVector.size(); ++i){
            showerToHitsMap[allShowerVector[i]] = hits_per_shower.at(allShowerVector[i].key());
        }

        art::FindOneP<recob::PFParticle> pfparticle_per_shower(allShowerHandle, evt, sss3dlabel);
        std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle> > showerToPFParticleMap;
        for(size_t i=0; i< allShowerVector.size(); ++i){
            showerToPFParticleMap[allShowerVector[i]] = pfparticle_per_shower.at(allShowerVector[i].key());
        }

        art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>("pandora");//""pandoraPatRec:allOutcomes");
        std::vector<art::Ptr<recob::PFParticle>> allPFParticleVector;
        art::fill_ptr_vector(allPFParticleVector,pfParticleHandle);

        //art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  "pandora");//PatRec:allOutcomes");
        //pfPartToMetadataAssoc.at(pfp.key());

        size_t n_all_shr = allShowerVector.size();
        m_sss3d_num_showers = (int)n_all_shr-showers.size();

        if(showers.size()==0) return;

        auto primary_shower = showers.front();

        std::cout<<"PandoraAllOutcomesShower has "<<n_all_shr<<" Showers in total. "<<std::endl;
        for(auto &shr: allShowerVector){
            //lets look at 3D distance to "vertex", and only go further with things that are within 80cm [default]
            double dist = sqrt(pow(m_vertex_pos_x - shr->ShowerStart().X(),2)+pow(m_vertex_pos_y - shr->ShowerStart().Y(),2)+pow(m_vertex_pos_z - shr->ShowerStart().Z(),2) );
            if(dist>max_conv_dist) continue;

            auto pfp = showerToPFParticleMap[shr];  
            //for(auto &prr: allPFParticleVector){
            //        std::cout<<pfp->Self()<<" "<<pfp.key()<<" "<<prr->Self()<<" "<<prr.key()<<std::endl;
            //}

            //What are we intested in learning
            std::vector<std::string> interested = {"IsClearCosmic","TrackScore","NuScore","IsNeutrino","SliceIndex"};

            /*
               std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatas = pfPartToMetadataAssoc.at(pfp.key());
               for(auto &meta: metadatas){
               std::map<std::string, float> propertiesmap  = meta->GetPropertiesMap();

               for (auto it:propertiesmap ){
               std::cout<<it.first<<" "<<it.second<<" ";
               }
               std::cout<<std::endl;

            //for each of the things in the list
            for(auto &s: interested){
            if(propertiesmap.count(s)==1){
            std::cout<<" "<<s<<" :  "<<propertiesmap[s]<<" ";
            }else{
            std::cout<<" NO "<<s<<" . ";
            }
            }
            }
            */

            //std::cout<<shr.key()<<" Length "<<shr->Length()<<" Dist "<<dist<<" Impact: "<<impact_paramater_shr(m_vertex_pos_x, m_vertex_pos_y, m_vertex_pos_z, shr)<<" pfp key: "<<showerToPFParticleMap[shr].key()<<" Self "<<showerToPFParticleMap[shr]->Self()<<std::endl;

            //OK we need to "remove" the  showers that are neutrino showers as well as those that are the "track"

            bool is_matched = false;

            for(auto &s: showers){
                //const art::Ptr<recob::Slice> s_slice = slice_per_pfparticle.at(NormalShowerToPFParticleMap[s].key());   
                //std::cout<<s.key()<<"shr is in slice "<<s_slice->ID()<<std::endl;
                if(s->ShowerStart().X() == shr->ShowerStart().X() || pfp->Self()== NormalShowerToPFParticleMap[s]->Self()){
                    //std::cout<<"Its a match!"<<std::endl;
                    is_matched = true;
                }
            }

            for(auto &s: tracks){
                //const art::Ptr<recob::Slice> s_slice = slice_per_pfparticle.at(NormalTrackToPFParticleMap[s].key());   
                //std::cout<<s.key()<<"trk is in slice "<<s_slice->ID()<<std::endl;
                if(pfp->Self()== NormalTrackToPFParticleMap[s]->Self()){
                    //std::cout<<"Its a match!"<<std::endl;
                    is_matched = true;
                }
            }

            if(is_matched) 
            {
                // std::cout<<"matched and continuing"<<std::endl; 
                continue;
            }

            double senergy = this->CalcEShower(showerToHitsMap[shr]); 
            double invar = implied_invar_mass(m_vertex_pos_x, m_vertex_pos_y, m_vertex_pos_z,  primary_shower, m_reco_shower_energy_max[0], shr, senergy);
            double implied_invar = invar_mass(primary_shower, m_reco_shower_energy_max[0], shr, senergy) ;
            double shr_score = 0.0; //need pfp and metadata to get score, and might give slice! (This will be harder..) but on reflection, kinda important. PCA spread might be a good rplacement.
            int is_clear_cosmic_slice = 0 ;
            int is_nu_slice = 0;


            m_sss3d_shower_start_x.push_back(shr->ShowerStart().X());
            m_sss3d_shower_start_y.push_back(shr->ShowerStart().Y());
            m_sss3d_shower_start_z.push_back(shr->ShowerStart().Z());
            m_sss3d_shower_dir_x.push_back(shr->Direction().X());
            m_sss3d_shower_dir_y.push_back(shr->Direction().Y());
            m_sss3d_shower_dir_z.push_back(shr->Direction().Z());
            m_sss3d_shower_length.push_back(shr->Length());
            m_sss3d_shower_conversion_dist.push_back(dist);
            m_sss3d_shower_invariant_mass.push_back(invar);
            m_sss3d_shower_implied_invariant_mass.push_back(implied_invar);
            double imp = impact_paramater_shr(m_vertex_pos_x, m_vertex_pos_y, m_vertex_pos_z, shr);
            m_sss3d_shower_impact_parameter.push_back(imp);

            if(dist!=0) {
                m_sss3d_shower_ioc_ratio.push_back(imp/dist);
            }else{
                m_sss3d_shower_ioc_ratio.push_back(0);

            }
            m_sss3d_shower_energy_max.push_back(senergy);// 
            m_sss3d_shower_score.push_back(shr_score);
            m_sss3d_slice_clear_cosmic.push_back(is_clear_cosmic_slice);
            m_sss3d_slice_nu.push_back(is_nu_slice);
        }   

        return;
    }



    void SinglePhoton::SimpleSecondShowerCluster(){

        std::string base = "sss3d_";
        std::vector<std::string> mod = {"ioc_ranked","invar_ranked"};

        m_sss3d_ioc_ranked_en = -9;
        m_sss3d_ioc_ranked_conv = -9;
        m_sss3d_ioc_ranked_invar = -9;
        m_sss3d_ioc_ranked_implied_invar = -9;
        m_sss3d_ioc_ranked_ioc = -9;
        m_sss3d_ioc_ranked_opang = -9;
        m_sss3d_ioc_ranked_implied_opang = -9; 
        m_sss3d_ioc_ranked_id = -9;

        m_sss3d_invar_ranked_en = -9;
        m_sss3d_invar_ranked_conv = -9;
        m_sss3d_invar_ranked_invar = -9;
        m_sss3d_invar_ranked_implied_invar = -9;
        m_sss3d_invar_ranked_ioc = -9;
        m_sss3d_invar_ranked_opang = -9;
        m_sss3d_invar_ranked_implied_opang = -9; 
        m_sss3d_invar_ranked_id = -9;


        std::string base2d = "sss_";
        std::vector<std::string> mod2d = {"ioc_ranked","conv_ranked","invar_ranked"};

        m_sss2d_ioc_ranked_en = -9;
        m_sss2d_ioc_ranked_conv = -9;
        m_sss2d_ioc_ranked_ioc = -9;
        m_sss2d_ioc_ranked_pca = -9;
        m_sss2d_ioc_ranked_invar = -9;
        m_sss2d_ioc_ranked_angle_to_shower = -9;
        m_sss2d_ioc_ranked_num_planes = -9;

        m_sss2d_conv_ranked_en = -9;
        m_sss2d_conv_ranked_conv = -9;
        m_sss2d_conv_ranked_ioc = -9;
        m_sss2d_conv_ranked_pca = -9;
        m_sss2d_conv_ranked_invar = -9;
        m_sss2d_conv_ranked_angle_to_shower = -9;
        m_sss2d_conv_ranked_num_planes = -9;

        m_sss2d_invar_ranked_en = -9;
        m_sss2d_invar_ranked_conv = -9;
        m_sss2d_invar_ranked_ioc = -9;
        m_sss2d_invar_ranked_pca = -9;
        m_sss2d_invar_ranked_invar = -9;
        m_sss2d_invar_ranked_angle_to_shower = -9;
        m_sss2d_invar_ranked_num_planes = -9;

        //---------------------------------------
        //First off, the 3D showers 
            //First some 3D shower information
            if(m_sss3d_shower_conversion_dist.size()>0 && m_reco_shower_energy_max.size()>0){
                //std::cout<<"Primary shower en "<<reco_shower_energy_max->at(0)<<std::endl;

                std::vector<double> inv = m_sss3d_shower_implied_invariant_mass;
                for(auto &v : inv) v = fabs(v-m_mass_pi0_mev);

                std::vector<size_t> ranked_ioc = sort_indexes_rev<double>((m_sss3d_shower_ioc_ratio));
                std::vector<size_t> ranked_invar = sort_indexes_rev<double>((inv));
                std::vector<size_t> ranked_conv = sort_indexes_rev<double>((m_sss3d_shower_conversion_dist));
                std::vector<size_t> ranked_en = sort_indexes_rev<double>((m_sss3d_shower_energy_max));

                int to_consider = m_sss3d_shower_conversion_dist.size();

                if(false){
                    std::cout<<"IOC"<<std::endl;
                    for(int j=0; j<to_consider; j++){
                        std::cout<<"--"<<ranked_ioc[j]<<" ioc: "<<m_sss3d_shower_ioc_ratio.at( ranked_ioc[j] )<<" invar: "<<m_sss3d_shower_implied_invariant_mass.at(ranked_ioc[j])<<" en: "<<m_sss3d_shower_energy_max.at(ranked_ioc[j])<<" conv: "<<m_sss3d_shower_conversion_dist.at(ranked_ioc[j])<<std::endl;
                    }

                    std::cout<<"INVAR"<<std::endl;
                    for(int j=0; j<to_consider; j++){
                        std::cout<<"--"<<ranked_invar[j]<<" ioc: "<<m_sss3d_shower_ioc_ratio.at( ranked_invar[j] )<<" invar: "<<m_sss3d_shower_implied_invariant_mass.at(ranked_invar[j])<<" en: "<<m_sss3d_shower_energy_max.at(ranked_invar[j])<<" conv: "<<m_sss3d_shower_conversion_dist.at(ranked_invar[j]) <<std::endl;
                    }

                    std::cout<<"EN"<<std::endl;
                    for(int j=0; j<to_consider; j++){
                        int rk = ranked_en[ranked_en.size()-1-j];
                        std::cout<<"--"<<rk<<" ioc: "<<m_sss3d_shower_ioc_ratio.at( rk )<<" invar: "<<m_sss3d_shower_implied_invariant_mass.at(rk)<<" en: "<<m_sss3d_shower_energy_max.at(rk)<<" conv: "<<m_sss3d_shower_conversion_dist.at(rk)<<std::endl;
                    }
                    std::cout<<"CONV"<<std::endl;
                    for(int j=0; j<to_consider; j++){
                        std::cout<<"--"<<ranked_conv[j]<<" ioc: "<<m_sss3d_shower_ioc_ratio.at( ranked_conv[j] )<<" invar: "<<m_sss3d_shower_implied_invariant_mass.at(ranked_conv[j])<<" en: "<<m_sss3d_shower_energy_max.at(ranked_conv[j])<<" conv: "<<m_sss3d_shower_conversion_dist.at(ranked_conv[j])<<std::endl;
                    }
                }


                //std::cout<<"Best IOC "<<"--"<<ranked_ioc[0]<<" ioc: "<<sss3d_shower_ioc_ratio->at( ranked_ioc[0] )<<" invar: "<<sss3d_shower_implied_invariant_mass->at(ranked_ioc[0])<<" en: "<<sss3d_shower_energy_max->at(ranked_ioc[0])<<" conv: "<<sss3d_shower_conversion_dist->at(ranked_ioc[0])<<std::endl;
                m_sss3d_ioc_ranked_en = m_sss3d_shower_energy_max.at(ranked_ioc[0]);            
                m_sss3d_ioc_ranked_conv = m_sss3d_shower_conversion_dist.at(ranked_ioc[0]);            
                m_sss3d_ioc_ranked_invar = m_sss3d_shower_invariant_mass.at(ranked_ioc[0]);            
                m_sss3d_ioc_ranked_implied_invar = m_sss3d_shower_implied_invariant_mass.at(ranked_ioc[0]);            
                m_sss3d_ioc_ranked_ioc = m_sss3d_shower_ioc_ratio.at(ranked_ioc[0]);
                m_sss3d_ioc_ranked_opang = 1.0 - pow(m_sss3d_shower_invariant_mass.at(ranked_ioc[0]),2)/(2.0*m_sss3d_shower_energy_max.at(ranked_ioc[0])*m_reco_shower_energy_max.at(0));
                m_sss3d_ioc_ranked_implied_opang = 1.0 - pow(m_sss3d_shower_implied_invariant_mass.at(ranked_ioc[0]),2)/(2.0*m_sss3d_shower_energy_max.at(ranked_ioc[0])*m_reco_shower_energy_max.at(0));
                m_sss3d_ioc_ranked_id = ranked_ioc[0];


                // std::cout<<"Best invar "<<"--"<<ranked_invar[0]<<" ioc: "<<sss3d_shower_ioc_ratio.at( ranked_invar[0] )<<" invar: "<<sss3d_shower_implied_invariant_mass.at(ranked_invar[0])<<" en: "<<sss3d_shower_energy_max.at(ranked_invar[0])<<" conv: "<<sss3d_shower_conversion_dist.at(ranked_invar[0])<<std::endl;

		// minimum discrepancy between implied invariant mass and pi0 mass
                m_sss3d_invar_ranked_en = m_sss3d_shower_energy_max.at(ranked_invar[0]);            
                m_sss3d_invar_ranked_conv = m_sss3d_shower_conversion_dist.at(ranked_invar[0]);            
                m_sss3d_invar_ranked_invar = m_sss3d_shower_invariant_mass.at(ranked_invar[0]);            
                m_sss3d_invar_ranked_implied_invar = m_sss3d_shower_implied_invariant_mass.at(ranked_invar[0]);            
                m_sss3d_invar_ranked_ioc = m_sss3d_shower_ioc_ratio.at(ranked_invar[0]);
                m_sss3d_invar_ranked_opang = 1.0 - pow(m_sss3d_shower_invariant_mass.at(ranked_invar[0]),2)/(2.0*m_sss3d_shower_energy_max.at(ranked_invar[0])*m_reco_shower_energy_max.at(0));
                m_sss3d_invar_ranked_implied_opang = 1.0 - pow(m_sss3d_shower_implied_invariant_mass.at(ranked_invar[0]),2)/(2.0*m_sss3d_shower_energy_max.at(ranked_invar[0])*m_reco_shower_energy_max.at(0));
                m_sss3d_invar_ranked_id = ranked_invar[0];

            }//end of 3D shower searching


            //Now some 2D shower information
            //
            if(m_sss_num_candidates>0){
                //std::cout<<"2D clusters: "<<sss_num_candidates<<std::endl;
                std::vector<int> nplans(3,0);
                std::vector<std::vector<int>> indexmap(3);


                for(int i=0; i< m_sss_num_candidates; i++){
                    //std::cout<<i<<" p: "<<m_sss_candidate_plane->at(i)<<" pdg: "<<m_sss_candidate_parent_pdg->at(i)<<" ovf "<<m_sss_candidate_overlay_fraction->at(i)<<" conv: "<<m_sss_candidate_min_dist->at(i)<<std::endl;
                }

                std::vector<std::vector<int>> uniq_candidates;

                for(int i=0; i< m_sss_num_candidates; i++){
                    int ip = m_sss_candidate_plane.at(i);
                    //int nhits = sss_candidate_num_hits.at(i);
                    nplans[ip]++;
                    indexmap[ip].push_back(i);

                    //Two passes to build up all "Candidates" for 2 and 3 plane matches
                    for(int j=i;j<m_sss_num_candidates;j++){
                        int jp = m_sss_candidate_plane.at(j);
                        if(jp==ip) continue;

                        bool contain_ij = false;
                        bool contain_ji = false;
                        if(m_sss_candidate_mean_tick.at(j)<=m_sss_candidate_max_tick.at(i) && m_sss_candidate_mean_tick.at(j) >= m_sss_candidate_min_tick.at(i))contain_ij = true;
                        if(m_sss_candidate_mean_tick.at(i)<=m_sss_candidate_max_tick.at(j) && m_sss_candidate_mean_tick.at(i) >= m_sss_candidate_min_tick.at(j))contain_ji = true;
                        //                                std::cout<<i<<" "<<j<<" "<<contain_ij<<" "<<contain_ji<<std::endl;
                        if(contain_ij && contain_ji){
                            uniq_candidates.push_back({i,j});
                        }
                    }
                }

                //Now loop over to check if any indlude a third plane
                for(int i = 0; i< (int)uniq_candidates.size(); i++){
                    for(int k=0; k<m_sss_num_candidates; k++){
                        //first check if this possible 3rd match is on a seperate plane
                        bool new_plane = true;
                        for(auto &pp:uniq_candidates[i]){
                            if(m_sss_candidate_plane.at(k)==m_sss_candidate_plane.at(pp) ) new_plane = false;   
                        }
                        if(new_plane){

                            bool contain_ik = false;
                            bool contain_ki = false;
                            bool contain_jk = false;
                            bool contain_kj = false;
                            if(m_sss_candidate_mean_tick.at(k)<=m_sss_candidate_max_tick.at(uniq_candidates[i][0]) && m_sss_candidate_mean_tick.at(k) >= m_sss_candidate_min_tick.at(uniq_candidates[i][0]))contain_ik = true;
                            if(m_sss_candidate_mean_tick.at(uniq_candidates[i][0])<=m_sss_candidate_max_tick.at(k) && m_sss_candidate_mean_tick.at(uniq_candidates[i][0]) >= m_sss_candidate_min_tick.at(k))contain_ki = true;
                            if(m_sss_candidate_mean_tick.at(k)<=m_sss_candidate_max_tick.at(uniq_candidates[i][1]) && m_sss_candidate_mean_tick.at(k) >= m_sss_candidate_min_tick.at(uniq_candidates[i][1]))contain_ik = true;
                            if(m_sss_candidate_mean_tick.at(uniq_candidates[i][1])<=m_sss_candidate_max_tick.at(k) && m_sss_candidate_mean_tick.at(uniq_candidates[i][1]) >= m_sss_candidate_min_tick.at(k))contain_ki = true;

                            //If this matches well with Either last candidate, include as a possibility
                            if((contain_ik&&contain_ki) || (contain_jk&&contain_kj)){
                                uniq_candidates[i].push_back(k);
                            }

                        }
                    }
                }
                //Check which candidates have been used where
                std::vector<int> used_candidates(m_sss_num_candidates);
                for(int i = 0; i< (int)uniq_candidates.size(); i++){
                    for(auto &j: uniq_candidates[i]){
                        used_candidates[j]++;
                    }
                }

                //If a candidate has been included in NO 2 or 3 plane cluster, treat it on its own
                for(int i = 0; i< (int)used_candidates.size(); i++){
                    if(used_candidates[i]==0) uniq_candidates.push_back({i});
                }

                //Now lets delete any permutations
                std::vector<std::vector<int>> uniq_candidates2;
                uniq_candidates2.push_back(uniq_candidates.front()); 

                for(int i = 1; i< (int)uniq_candidates.size(); i++){

                    bool perm = false;
                    for(int j = 0; j< (int)uniq_candidates2.size(); j++){
                        perm = marks_compare_vec_nonsense<int>(uniq_candidates[i], uniq_candidates2[j]);
                        if(perm) break;
                    }
                    if(!perm) uniq_candidates2.push_back(uniq_candidates[i]);
                }

                //Printing candidates (After perm check)
                std::cout<<"After: used_candidates "<<m_sss_num_candidates<<std::endl;
                for(int i = 0; i< (int)uniq_candidates2.size(); i++){
                    std::cout<<i<<" | ";
                    for(auto &j: uniq_candidates2[i])std::cout<<" "<<j;
                    std::cout<<std::endl;
                }

                //Right lets CULL and rank the candidates
                std::vector<bool> candidate_pass(uniq_candidates2.size(),false);
                std::vector<double>  candidates_en(uniq_candidates2.size(),0);
                std::vector<double>  candidates_ioc(uniq_candidates2.size(),0);
                std::vector<double>  candidates_conv(uniq_candidates2.size(),0);
                std::vector<double>  candidates_pca(uniq_candidates2.size(),0);
                std::vector<double>  candidates_angle_to_shower(uniq_candidates2.size(),0);
                std::vector<int>     candidates_num_planes(uniq_candidates2.size(),0);
                std::vector<double>  candidates_eff_invar(uniq_candidates2.size(),0);
                std::vector<double>  candidates_eff_invar_diff(uniq_candidates2.size(),0);

                //rank by min_impat/max_min_dist and select
                //rank by Energy energy

                for(int j=0; j<(int)uniq_candidates2.size();j++){
                    int nt=uniq_candidates2[j].size();
                    //std::cout<<"Candidate #: "<<j<<" has "<<nt<<" 2D clusters"<<std::endl;

                    double mean_min_dist = 0.0;
                    double max_min_dist = 0.0;

                    double mean_energy = 0.0;

                    double mean_impact = 0.0;
                    double mean_conv = 0.0;
                    double min_conv = 999;

                    double min_impact = 999;
                    double mean_invar = 0.0;
                    double mean_invar_diff = 0.0;

                    double max_pca = 0;
                    double min_angle = 999;

        
                    //int is_min_slice = 999;
                    //std::vector<int> is_in_slice; 

                    for(int c=0; c< nt;++c){
                        int ic = uniq_candidates2[j][c];

                        //std::cout<<"----- plane: "<<m_sss_candidate_plane.at(ic)<<" nhits: "<<m_sss_candidate_num_hits.at(ic)<<" imp: "<<m_sss_candidate_impact_parameter.at(ic)<<" en: "<<m_sss_candidate_energy.at(ic)<<" a2s: "<<m_sss_candidate_angle_to_shower.at(ic)<<" conv: "<<m_sss_candidate_min_dist.at(ic)<<" pdg: "<<m_sss_candidate_parent_pdg.at(ic)<<std::endl;

                        double eff_invar = sqrt(2.0*m_sss_candidate_energy.at(ic)*m_reco_shower_energy_max.at(0)*(1.0-cos(m_sss_candidate_angle_to_shower.at(ic))));
                        double eff_invar_diff = fabs(eff_invar-139.5);

                        //is_min_slice = std::min(is_in_slice, m_sss_candidate_in_nu_slice[ic] );
                        //is_in_slice.push_back(m_sss_candidate_in_nu_slice[ic]); 

                        mean_min_dist +=m_sss_candidate_min_dist.at(ic)/(double)nt;
                        mean_energy +=m_sss_candidate_energy.at(ic)/(double)nt;
                        mean_impact +=m_sss_candidate_impact_parameter.at(ic)/(double)nt;
                        mean_conv +=m_sss_candidate_min_dist.at(ic)/(double)nt;
                        mean_invar +=eff_invar/(double)nt;
                        mean_invar_diff +=eff_invar_diff/(double)nt;

                        min_conv = std::min(min_conv, m_sss_candidate_min_dist.at(ic));
                        max_min_dist = std::max(max_min_dist, m_sss_candidate_min_dist.at(ic));
                        max_pca = std::max(max_pca, m_sss_candidate_PCA.at(ic));
                        min_impact = std::min(min_impact, m_sss_candidate_impact_parameter.at(ic));
                        min_angle = std::min(min_angle, m_sss_candidate_angle_to_shower.at(ic));

                    }
                    std::cout<<"======== Mean En "<<mean_energy<<" mean dist "<<mean_min_dist<<" meanimpact "<<mean_impact<<" minimpact/maxdist :  "<<min_impact/max_min_dist<<" invar "<<mean_invar<<std::endl;
                    candidates_ioc[j]=min_impact/max_min_dist;
                    candidates_en[j]=mean_energy;
                    candidates_conv[j] = min_conv;
                    candidates_pca[j] = max_pca;
                    candidates_angle_to_shower[j] = min_angle;
                    candidates_num_planes[j] =nt; 
                    candidates_eff_invar_diff[j] = mean_invar_diff;
                    candidates_eff_invar[j] = mean_invar;
                    //candidates_min_slice[j] = is_min_slice;
                    //candidates_in_slice[j] = is_in_slice;
                }

                std::vector<size_t> ranked_ioc = sort_indexes_rev<double>(candidates_ioc);
                std::vector<size_t> ranked_invar = sort_indexes_rev<double>(candidates_eff_invar_diff);
                std::vector<size_t> ranked_conv = sort_indexes_rev<double>(candidates_conv);

                std::cout<<"========== Ranking ======== "<<std::endl;
                std::cout<<"IOC ";   for (auto &ii: ranked_ioc) std::cout<<" "<<ii; std::cout<<std::endl;
                std::cout<<"CONV ";   for (auto &ii: ranked_conv) std::cout<<" "<<ii; std::cout<<std::endl;
                std::cout<<"INVAR ";for (auto &ii: ranked_invar) std::cout<<" "<<ii; std::cout<<std::endl;

                m_sss2d_ioc_ranked_en = candidates_en[ranked_ioc[0]];
                m_sss2d_ioc_ranked_conv = candidates_conv[ranked_ioc[0]];
                m_sss2d_ioc_ranked_ioc = candidates_ioc[ranked_ioc[0]];
                m_sss2d_ioc_ranked_invar = candidates_eff_invar[ranked_ioc[0]];
                m_sss2d_ioc_ranked_pca = candidates_pca[ranked_ioc[0]];
                m_sss2d_ioc_ranked_angle_to_shower  = candidates_angle_to_shower[ranked_ioc[0]];
                m_sss2d_ioc_ranked_num_planes  = candidates_num_planes[ranked_ioc[0]];

                m_sss2d_conv_ranked_en = candidates_en[ranked_conv[0]];
                m_sss2d_conv_ranked_conv = candidates_conv[ranked_conv[0]];
                m_sss2d_conv_ranked_ioc = candidates_ioc[ranked_conv[0]];
                m_sss2d_conv_ranked_invar = candidates_eff_invar[ranked_conv[0]];
                m_sss2d_conv_ranked_pca = candidates_pca[ranked_conv[0]];
                m_sss2d_conv_ranked_angle_to_shower  = candidates_angle_to_shower[ranked_conv[0]];
                m_sss2d_conv_ranked_num_planes  = candidates_num_planes[ranked_conv[0]];

                m_sss2d_invar_ranked_en = candidates_en[ranked_invar[0]];
                m_sss2d_invar_ranked_conv = candidates_conv[ranked_invar[0]];
                m_sss2d_invar_ranked_ioc = candidates_ioc[ranked_invar[0]];
                m_sss2d_invar_ranked_invar = candidates_eff_invar[ranked_invar[0]];
                m_sss2d_invar_ranked_pca = candidates_pca[ranked_invar[0]];
                m_sss2d_invar_ranked_angle_to_shower  = candidates_angle_to_shower[ranked_invar[0]];
                m_sss2d_invar_ranked_num_planes  = candidates_num_planes[ranked_invar[0]];
            }

        return ;
    }


   
    std::pair<bool, std::vector<double>> SinglePhoton::clusterCandidateOverlap(const std::vector<int> & candidate_indices, const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks){

        size_t size = candidate_indices.size();
	if(size == 0){
	    throw std::runtime_error("SinglePhoton::clusterCandidateOverlap: No cluster candidates to analyze time overlap for..");
	}

	// at most 3 cluster indices (for 3 planes)
	std::vector<int> planes;
        std::vector<double> max_ticks;
        std::vector<double> min_ticks;
        std::vector<double> tick_length;

        for(auto i : candidate_indices){
            planes.push_back(cluster_planes[i]);
            
            max_ticks.push_back(cluster_max_ticks[i]);
            min_ticks.push_back(cluster_min_ticks[i]); 
            tick_length.push_back(cluster_max_ticks[i] - cluster_min_ticks[i]);
        }


	//if candidates are not on different planes
	if( size == 2 && planes[0] == planes[1])
            return {false, std::vector<double>(2, -1.0)};
        if( size == 3 && (planes[0] == planes[1] || planes[1] == planes[2] || planes[0] == planes[2]))
            return {false, std::vector<double>(3, -1.0)};

	//calculate the overlapping tick-span
	double tick_overlap = DBL_MAX;

	//can be simplied as picking the minimum max_tick and maximum min_tick and do the subtraction
        for(auto max_e : max_ticks)
            for(auto min_e : min_ticks)
                if(max_e - min_e < tick_overlap)
                    tick_overlap = max_e - min_e;
 
	// if tick overlap is negative, meaning these clusters are not overlapping
   	if(tick_overlap < 0)
            return {false, std::vector<double>(size, -1.0)};
        else{
            std::vector<double> overlap_fraction;
            for(auto l: tick_length){
                overlap_fraction.push_back( l==0? 1.0 : tick_overlap/l);
	    }
            return {true, overlap_fraction};
        }
    }

   
    std::pair<int, std::pair<std::vector<std::vector<double>>, std::vector<double>>> SinglePhoton::GroupClusterCandidate(int num_clusters,  const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks){
	std::cout << "SinglePhoton::group_cluster_candidate\t|| Total of " << num_clusters << " to be grouped" << std::endl;

	int num_cluster_groups=0; // number of matched cluster groups in total
	std::vector<std::vector<double>> grouped_cluster_indices;
	std::vector<double> cluster_group_timeoverlap_fraction;
	if(num_clusters <= 1)
	    return {num_cluster_groups, {grouped_cluster_indices, cluster_group_timeoverlap_fraction}};

        for(int i = 0; i != num_clusters -1; ++i){
            for(int j = i+1; j != num_clusters; ++j){

		//first, look at candidate pairs
                auto pair_result = clusterCandidateOverlap({i,j}, cluster_planes, cluster_max_ticks, cluster_min_ticks);
                if( pair_result.first){

		    ++num_cluster_groups;
		    grouped_cluster_indices.push_back({(double)i,(double)j});
                    double min_frac = *std::min_element(pair_result.second.cbegin(), pair_result.second.cend());
		    cluster_group_timeoverlap_fraction.push_back(min_frac);
		    std::cout << "Grouped cluster candidate: (" << i  << ", " << j << ") | Minimum time tick overlap fraction: " << min_frac << std::endl;

		    // if the pair is succefully grouped, look at possible trios
                    for(int k = j+1; k!= num_clusters; ++k){
                        auto tri_result = clusterCandidateOverlap({i,j,k}, cluster_planes, cluster_max_ticks, cluster_min_ticks);
                        if(tri_result.first){
			    ++num_cluster_groups;
                    	    grouped_cluster_indices.push_back({(double)i,(double)j,(double)k});
                            min_frac = *std::min_element(tri_result.second.cbegin(), tri_result.second.cend());
			    cluster_group_timeoverlap_fraction.push_back(min_frac);
			    std::cout << "Grouped cluster candidate: (" << i  << ", " << j << ", " << k << ") | Minimum time tick overlap fraction: " << min_frac << std::endl;
                        }
                    } //k loop
                }
            }//j loop
        }//i loop

	std::cout << "SinglePhoton::GroupClusterCandidate\t|| Formed " << num_cluster_groups << " cluster groups" << std::endl;

	return {num_cluster_groups, {grouped_cluster_indices, cluster_group_timeoverlap_fraction}};
    } 

}
