#include "SinglePhoton_module.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TEllipse.h"
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



    void SinglePhoton::SecondShowerSearch(
            const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
            const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
            const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
            const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
            std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap,
					detinfo::DetectorPropertiesData const & theDetector
			){



        std::cout<<"ERROR! SecondShowerSearch has been made redundant by SEAview methodology. see SEAview/SEAview.h for details. This shouldnt be running at all"<<std::endl;
        std::cerr<<"ERROR! SecondShowerSearch has been made redundant by SEAview methodology. see SEAview/SEAview.h for details. This shouldnt be running at all"<<std::endl;
        exit(EXIT_FAILURE);

        int total_track_hits =0;
        int total_shower_hits =0;
        int nu_slice_id = -999;

        std::vector<art::Ptr<recob::Hit>> associated_hits;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits_plane0;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits_plane1;
        std::vector<art::Ptr<recob::Hit>> unassociated_hits_plane2;

        std::vector< std::map<size_t, std::vector<art::Ptr<recob::Hit>>>> v_newClusterToHitsMap(3);//one for each plane

        for(size_t t =0; t< tracks.size(); t++){
            art::Ptr<recob::Track> track = tracks[t];
            art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];
            int sliceid = pfParticleToSliceIDMap.at(pfp);
            auto slicehits = sliceIDToHitsMap.at(sliceid);
            auto trackhits = pfParticleToHitsMap.at(pfp);

            std::cout<<"SinglePhoton::SSS\t||\ttrack "<<t<<" is in slice "<<sliceid<<" which has "<<slicehits.size()<<" hits. This track has  "<<trackhits.size()<<" of them. "<<std::endl;
            total_track_hits+=trackhits.size();
            if(nu_slice_id !=  sliceid && nu_slice_id != -999){
                std::cout<<"ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
                exit(EXIT_FAILURE);
            }   
            nu_slice_id = sliceid;


            for(auto &h: trackhits){
                associated_hits.push_back(h);
            }
        }

        for(size_t s =0; s< showers.size(); s++){
            art::Ptr<recob::Shower> shower = showers[s];
            art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(shower);
            int sliceid = pfParticleToSliceIDMap.at(pfp);
            auto slicehits = sliceIDToHitsMap.at(sliceid);
            auto showerhits = pfParticleToHitsMap.at(pfp);

            std::cout<<"SinglePhoton::SSS\t||\tshower "<<s<<" is in slice "<<sliceid<<" which has "<<slicehits.size()<<" hits. This shower has  "<<showerhits.size()<<" of them. "<<std::endl;
            total_shower_hits+=showerhits.size();
            if(nu_slice_id !=  sliceid && nu_slice_id!=-999){
                std::cout<<"ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
                exit(EXIT_FAILURE);
            }
            nu_slice_id = sliceid;

            for(auto &h: showerhits){
                associated_hits.push_back(h);
            }

        }

        std::cout<<"SinglePhoton::SSS\t||\tSo in total we have "<<total_shower_hits<<" shower hits and "<<total_track_hits<<" track hits"<<" assocaedHits total: "<<associated_hits.size()<<std::endl;
        m_sss_num_associated_hits = total_shower_hits+total_track_hits;

        if(nu_slice_id>=0){
            std::cout<<"SinglePhoton::SSS\t||\t So that leaves "<<sliceIDToHitsMap.at(nu_slice_id).size()-total_shower_hits-total_track_hits<<" hits not included in tracks and showers"<<std::endl;
            m_sss_num_unassociated_hits = sliceIDToHitsMap.at(nu_slice_id).size()-total_shower_hits-total_track_hits;

            if(m_sss_num_unassociated_hits<0){
                std::cout<<"ERROR!! Number of unassociated hits is negative, i.e: num_associated: "<<m_sss_num_associated_hits<<" and total slice hits: "<<sliceIDToHitsMap.at(nu_slice_id).size()<<std::endl;
                exit(EXIT_FAILURE);
            }


            std::vector<art::Ptr<recob::Hit>> slicehits = sliceIDToHitsMap.at(nu_slice_id);
            for(auto &h: slicehits){

                bool is_associated = false;
                for(auto &a: associated_hits){
                    if(h==a){
                        is_associated = true;
                        break;
                    }
                }

                if(!is_associated){
                    unassociated_hits.push_back(h);
                    auto plane_view = h->View();
                    switch((int)plane_view){
                        case (0) :
                            unassociated_hits_plane0.push_back(h);
                            break;
                        case (1) :
                            unassociated_hits_plane1.push_back(h);
                            break;
                        case (2) :
                            unassociated_hits_plane2.push_back(h);
                            break;
                    }

                }

            }

            std::vector<std::vector<art::Ptr<recob::Hit>>> unassociated_hits_all = {unassociated_hits_plane0,unassociated_hits_plane1,unassociated_hits_plane2};

            std::cout<<"SinglePhoton::SSS\t||\tassociated_hits.size() "<<associated_hits.size()<<" unassociated_hits.size() "<<unassociated_hits.size()<<" p0: "<<unassociated_hits_plane0.size()<<" p1:  "<<unassociated_hits_plane1.size()<<" p2: "<<unassociated_hits_plane2.size()<<std::endl;


            if(bool_make_sss_plots && showers.size()==1 && tracks.size()>0){

                //TFile *f = new TFile("t.root","recreate");
                //f->cd();

                std::string print_name = "sss_"+std::to_string(m_run_number)+"_"+std::to_string(m_subrun_number)+"_"+std::to_string(m_event_number);
                TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,2400);
                can->Divide(4,3,0,0.1);

                double tick_max = 0;
                double tick_min = 1e10;
                std::vector<double> chan_max(3,0);
                std::vector<double> chan_min(3,1e10);


                //First grab all shower clusters
                std::vector<std::vector<TGraph *>> pts_shr( showers.size(), std::vector<TGraph *>(3)  );

                for(size_t s =0; s< showers.size(); s++){
                    art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(showers[s]);
                    auto showerhits = pfParticleToHitsMap.at(pfp);

                    std::vector<TGraph*> t_pts(3);
                    std::vector<std::vector<double>>  vec_t(3);
                    std::vector<std::vector<double>>  vec_c(3);

                    for(auto &h: showerhits){
                        double wire = (double)h->WireID().Wire;
                        vec_c[(int)h->View()].push_back(wire);
                        //vec_c[(int)h->View()].push_back((double)h->Channel());
                        vec_t[(int)h->View()].push_back((double)h->PeakTime());
                        tick_max = std::max(tick_max, (double)h->PeakTime());
                        tick_min = std::min(tick_min, (double)h->PeakTime());
                        chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                        chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);
                    }
                    t_pts[0] = new TGraph(vec_c[0].size(), &(vec_c[0])[0], &(vec_t[0])[0]);
                    t_pts[1] = new TGraph(vec_c[1].size(), &(vec_c[1])[0], &(vec_t[1])[0]);
                    t_pts[2] = new TGraph(vec_c[2].size(), &(vec_c[2])[0], &(vec_t[2])[0]);

                    pts_shr[s] = t_pts;
                }


                std::vector<std::vector<TGraph *>> pts_trk( tracks.size(), std::vector<TGraph *>(3)  );
                for(size_t t =0; t< tracks.size(); t++){
                    art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap.at(tracks[t]);
                    auto trackhits = pfParticleToHitsMap.at(pfp);

                    std::vector<TGraph*> t_pts(3);
                    std::vector<std::vector<double>>  vec_t(3);
                    std::vector<std::vector<double>>  vec_c(3);

                    for(auto &h: trackhits){
                        double wire = (double)h->WireID().Wire;
                        vec_c[(int)h->View()].push_back(wire);
                        vec_t[(int)h->View()].push_back((double)h->PeakTime());
                        tick_max = std::max(tick_max, (double)h->PeakTime());
                        tick_min = std::min(tick_min, (double)h->PeakTime());
                        chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                        chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

                    }
                    t_pts[0] = new TGraph(vec_c[0].size(), &(vec_c[0])[0], &(vec_t[0])[0]);
                    t_pts[1] = new TGraph(vec_c[1].size(), &(vec_c[1])[0], &(vec_t[1])[0]);
                    t_pts[2] = new TGraph(vec_c[2].size(), &(vec_c[2])[0], &(vec_t[2])[0]);

                    pts_trk[t] = t_pts;
                }
                //Now the "Unassociated Hits"

                std::vector<TGraph *> g_unass(3);
                std::vector<std::vector<std::vector<double>>> pts_to_recluster(3); //plane, point, {wire,tic}
            //make a map tp actual hits here I guess.
            std::vector<std::map<int, art::Ptr<recob::Hit>>> mapPointIndexToHit(3);

            for(int i=0; i<3; i++){

                std::vector<double> vec_t;
                std::vector<double> vec_c;

                for(auto &h: unassociated_hits_all[i]){

                    double wire = (double)h->WireID().Wire;
                    vec_c.push_back(wire);
                    vec_t.push_back((double)h->PeakTime());

                    tick_max = std::max(tick_max, (double)h->PeakTime());
                    tick_min = std::min(tick_min, (double)h->PeakTime());
                    chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
                    chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

                    //for reclustering
                    std::vector<double> pt = {wire,vec_t.back()};
                    pts_to_recluster[(int)h->View()].push_back(pt);
                    mapPointIndexToHit[(int)h->View()][pts_to_recluster[(int)h->View()].size()-1] = h;
                }

                g_unass[i] = new TGraph(vec_c.size(), &vec_c[0], &vec_t[0]);


            }
            //Plotting now
            double plot_point_size=0.6;        

            std::vector<int> tcols = {kRed-7, kBlue-7, kGreen-3, kOrange-3, kCyan-3, kMagenta-3, kGreen+1 , kRed+1};
            int used_col=0;
            if(showers.size()+tracks.size() > tcols.size()){
                for(int i =0; i< (int)(showers.size()+tracks.size() - tcols.size() +2); i++){
                    tcols.push_back(tcols[(int)rangen->Uniform(0,7)]+(int)rangen->Uniform(-5,5));
                }
            }

            std::cout<<"SinglePhoton::SSS\t||\tTick Min: "<<tick_min<<" Max: "<<tick_max<<std::endl;
            auto const TPC = (*geom).begin_TPC();
            auto ID = TPC.ID();
            int fCryostat = ID.Cryostat;
            int fTPC = ID.TPC;
            std::cout<<"SinglePhoton::SSS\t||\t" << TPC.ID() << "= the beginning TPC ID" <<std::endl;
            std::cout<<"SinglePhoton::SSS\t||\tthe cryostat id = " << fCryostat << std::endl;  
            std::cout<<"SinglePhoton::SSS\t||\tthe tpc id = " << fTPC << std::endl;  


            //Plotting the vertex position on the plot.

            std::vector<double> vertex_time(3); 
            std::vector<double> vertex_wire(3); 


            std::vector<TGraph*> g_vertex(3);
            for(int i=0; i<3; i++){
                TPad * pader = (TPad*)can->cd(i+1);

                if(i==0 || i ==4 || i == 8) pader->SetLeftMargin(0.1);

                std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
                std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, theDetector)};

                vertex_time[i] = time[0];
                vertex_wire[i] = wire[0];

                if(i==0) m_vertex_pos_wire_p0 = wire[0];
                if(i==1) m_vertex_pos_wire_p1 = wire[0];
                if(i==2) m_vertex_pos_wire_p2 = wire[0];
                m_vertex_pos_tick = time[0];

                chan_max[i] = std::max( chan_max[i],wire[0]);
                chan_min[i] = std::min( chan_min[i],wire[0]);

                g_vertex[i] = new TGraph(1,&wire[0],&time[0]);
                g_vertex[i]->SetMarkerStyle(29);
                g_vertex[i]->SetMarkerSize(4);
                g_vertex[i]->SetMarkerColor(kMagenta-3);
                g_vertex[i]->GetYaxis()->SetRangeUser(tick_min*0.98,tick_max*1.02);
                g_vertex[i]->GetXaxis()->SetLimits(chan_min[i]*0.98,chan_max[i]*1.02);
                g_vertex[i]->SetTitle(("Plane " +std::to_string(i)).c_str());
                g_vertex[i]->GetYaxis()->SetTitle("Peak Hit Time Tick");
                g_vertex[i]->GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
                g_vertex[i]->Draw("ap");

                if(i>0){
                    g_vertex[i]->GetYaxis()->SetLabelOffset(999);
                    g_vertex[i]->GetYaxis()->SetLabelSize(0);
                }

                can->cd(i+5);
                g_vertex[i]->Draw("ap");

                can->cd(i+9);
                g_vertex[i]->Draw("ap");


            }

            //******************************** DeadWireRegions********************************************
            for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                int badchan = bad_channel_list_fixed_mcc9[i].first;                                       
                int ok = bad_channel_list_fixed_mcc9[i].second;       

                if(ok>1)continue;
                auto hs = geom->ChannelToWire(badchan);

                //std::cout<<"KNK: "<<bc<<" "<<hs[0]<<" "<<result.start().X()<<" "<<result.start().Y()<<" "<<result.start().Z()<<" "<<result.end().X()<<" "<<result.end().Y()<<" "<<result.end().Z()<<std::endl; 
                int thisp = (int)hs[0].Plane;
                double bc = hs[0].Wire;

                //                        std::cout<<"WIRE "<<thisp<<" "<<bc<<" "<<hs.size()<<std::endl;

                if(chan_min[thisp]*0.98 < bc && bc < chan_max[thisp]*1.02 ){
                    can->cd(thisp+1);
                    TLine *l = new TLine(bc,tick_min*0.98,bc,tick_max*1.02);
                    l->SetLineColor(kGray+1);
                    l->Draw("same");
                    can->cd(thisp+5);
                    l->Draw("same");
                    can->cd(thisp+9);
                    l->Draw("same");
                }
            }

            for(size_t t=0; t< pts_trk.size(); t++){
                int tcol = tcols[used_col];
                used_col++;

                for(int i=0; i<3; i++){
                    can->cd(i+1);
                    if(pts_trk[t][i]->GetN()>0){//need a check in case this track has no hits on this plane.
                        pts_trk[t][i]->Draw("p same"); 
                        pts_trk[t][i]->SetMarkerColor(tcol);
                        pts_trk[t][i]->SetFillColor(tcol);
                        pts_trk[t][i]->SetMarkerStyle(20);
                        pts_trk[t][i]->SetMarkerSize(plot_point_size);
                    }
                }
            }

            for(size_t t=0; t< pts_shr.size(); t++){
                int tcol = tcols[used_col];
                used_col++;

                for(int i=0; i<3; i++){
                    can->cd(i+1);
                    if(pts_shr[t][i]->GetN()>0){
                        pts_shr[t][i]->Draw("p same"); //used in the vertex
                        pts_shr[t][i]->SetMarkerColor(tcol);
                        pts_shr[t][i]->SetFillColor(tcol);
                        pts_shr[t][i]->SetMarkerStyle(20);
                        pts_shr[t][i]->SetMarkerSize(plot_point_size);
                    }
                }
            }


            for(int i=0; i<3; i++){
                can->cd(i+1);

                if(g_unass[i]->GetN()>0){
                    g_unass[i]->Draw("p same");
                    g_unass[i]->SetMarkerColor(kBlack);
                    g_unass[i]->SetMarkerStyle(20);
                    g_unass[i]->SetMarkerSize(plot_point_size);
                }
                g_vertex[i]->Draw("p same");
                can->cd(i+5);
                g_vertex[i]->Draw("p same");
                can->cd(i+9);
                g_vertex[i]->Draw("p same");

                double rad_cm = 12.0;
                TEllipse * ell_p = new TEllipse(vertex_wire[i],vertex_time[i],rad_cm/0.3,rad_cm*25);
                ell_p->SetLineColor(kRed);
                ell_p->SetFillStyle(0);
                ell_p->Draw("same");

            }






            //*****************************DBSCAN***********************************
			//CHECK
//            int min_pts = m_SEAviewDbscanMinPts;
//            double eps = m_SEAviewDbscanEps;
//            std::vector<int> num_clusters(3,0);
//
//            std::vector<std::vector<TGraph*>> g_clusters(3);
//            std::vector<std::vector<int>> cluster_labels(3);
//
//
//            std::vector<cluster> vec_clusters;
//
//            for(int i=0; i<3; i++){
//
//                std::cout<<"SinglePhoton::SSS\t||\tStarting to run DBSCAN for plane: "<<i<<" has "<<pts_to_recluster[i].size()<<" pts to do using eps: "<<eps<<" and min_pts: "<<min_pts<<std::endl; 
//                DBSCAN ReCluster(eps,min_pts);
//                cluster_labels[i] =  ReCluster.Scan2D(pts_to_recluster[i]);
//
//                for(auto &c: cluster_labels[i]){
//                    num_clusters[i] = std::max(c,num_clusters[i]);
//                }
//                std::cout << "SinglePhoton::SSS\t||\tOn this plane, "<<m_event_number<<" DBSCAN found: "<<num_clusters[i]<<" clusters"<<std::endl;
//            }
//
//            //Now fill the clusters
//
//            for(int i=0; i<3; i++){
//
//                for(int c=0; c<num_clusters[i]+1; c++){
//
//                    std::vector<std::vector<double>> pts;
//                    std::vector<art::Ptr<recob::Hit>> hitz;
//                    for(size_t p=0; p< pts_to_recluster[i].size(); p++){
//                        if(cluster_labels[i][p] == 0) continue;//noise 
//                        if(cluster_labels[i][p] == c){
//
//                            pts.push_back(pts_to_recluster[i][p]);
//                            hitz.push_back(mapPointIndexToHit[i].at(p));
//                        }
//
//                    }
//                    if(hitz.size()!=0){
//                        vec_clusters.emplace_back(c,i,pts,hitz);
//                        std::cout<<"SinglePhoton::SSS\t||\t Cluster "<<c<<" has "<<hitz.size()<<" hitz on plane and "<<pts.size()<<"points "<<i<<std::endl;
//                    }
//
//                }
//            }
//
//
//
//            //Create final maps;
//            for(size_t i=0; i<3; i++){
//                for(size_t p=0; p<pts_to_recluster[i].size(); p++){
//
//                    art::Ptr<recob::Hit> h = mapPointIndexToHit[i].at(p);// Get the hit
//                    size_t cluster_label = cluster_labels[i][p];//get the cluster index, 0 = noise                                                
//
//                    //std::cout<<i<<" "<<p<<" "<<cluster_label<<" "<<v_newClusterToHitsMap[i].count(cluster_label)<<std::endl;
//
//                    if(v_newClusterToHitsMap[i].count(cluster_label)<1){
//                        std::vector<art::Ptr<recob::Hit>> t_hs= {h};
//                        v_newClusterToHitsMap[i][cluster_label] = t_hs;
//                    }else{
//                        v_newClusterToHitsMap[i].at(cluster_label).push_back(h);//add it to list
//                    }
//                }
//            }
//
//
//            //for(size_t i=0; i<3; i++){
//            //  for(int c=0; c<num_clusters[i]+1; c++){
//            //auto v_newClusterToHitsMap[i][c];
//            // }
//            // }
//
//
//            int max_n_clusters = std::max(num_clusters[0]+1, std::max(num_clusters[1]+1,num_clusters[2]+1))+2;
//            std::vector<int> cluster_colors(max_n_clusters,0);
//            std::vector<int> base_col = {632,416, 600, 400, 616,  432,  800, 820,  840, 860,  880, 900};
//
//            for(int j=0; j< max_n_clusters; j++){
//                int b = (int)rangen->Uniform(0,11);
//                int mod = (int)rangen->Uniform(-10,+3);
//
//                cluster_colors[j] = base_col[b]+mod;
//            }
//            int c_offset = 0;
//
//            //Step next, loop over and make plots again
//            for(int i=0; i<3; i++){
//                std::vector<std::vector<double>> vec_time(num_clusters[i]+1);
//                std::vector<std::vector<double>> vec_wire(num_clusters[i]+1);
//                std::vector<TGraph*> tmp_g_clusters(num_clusters[i]+1);
//
//                if(cluster_labels[i].size() != pts_to_recluster[i].size()){
//                    std::cout<<"SinglePhoton::SSS\t||\tERROR!! someting amiss cluster labels of size "<<cluster_labels[i].size()<<" and pts  in this plane "<<pts_to_recluster[i].size()<<std::endl;
//                }  
//
//                for(size_t k=0; k< pts_to_recluster[i].size(); k++){
//                    //std::cout<<vec_wire.size()<<" "<<cluster_labels[i][k]<<std::endl;
//                    vec_wire[cluster_labels[i][k]].push_back(pts_to_recluster[i][k][0]); 
//                    vec_time[cluster_labels[i][k]].push_back(pts_to_recluster[i][k][1]); 
//
//                }
//
//                for(int c=0; c< num_clusters[i]+1; c++){
//                    int tcol = kBlack;
//                    if(c>0) tcol = cluster_colors[c+c_offset];
//                    tmp_g_clusters[c] = new TGraph(vec_wire[c].size(),&(vec_wire[c])[0],&(vec_time[c])[0] );
//                    can->cd(i+5);
//                    if(
//                            tmp_g_clusters[c]->GetN()>0){
//                        tmp_g_clusters[c]->Draw("p same");
//                        tmp_g_clusters[c]->SetMarkerColor(tcol);
//                        tmp_g_clusters[c]->SetFillColor(tcol);
//                        tmp_g_clusters[c]->SetMarkerStyle(20);
//                        tmp_g_clusters[c]->SetMarkerSize(plot_point_size);
//                    }
//                }
//                g_clusters[i] = tmp_g_clusters;
//                c_offset += num_clusters[i];
//            }

            //******************* INFO Plotting *******************************

            TPad *p_top_info = (TPad*)can->cd(4);
            p_top_info->cd();

            TLatex pottex;
            pottex.SetTextSize(0.045);
            pottex.SetTextAlign(13);  //align at top
            pottex.SetNDC();
            std::string pot_draw = "Run: "+std::to_string(m_run_number)+" SubRun: "+std::to_string(m_subrun_number)+" Event: "+std::to_string(m_event_number);
            pottex.DrawLatex(.1,.94, pot_draw.c_str());

            TLegend * l_top = new TLegend(0.5,0.5,0.85,0.85);

            for(size_t t=0; t< pts_shr.size(); t++){
                std::string sname = "Shower "+std::to_string(t);
                if(pts_shr[t][0]->GetN()>0){
                    l_top->AddEntry(pts_shr[t][0],sname.c_str(),"f");
                }else if(pts_shr[t][1]->GetN()>0){
                    l_top->AddEntry(pts_shr[t][1],sname.c_str(),"f");
                }else if(pts_shr[t][2]->GetN()>0){
                    l_top->AddEntry(pts_shr[t][2],sname.c_str(),"f");
                }
            }

            for(size_t t=0; t< pts_trk.size(); t++){
                std::string sname = "Track "+std::to_string(t);
                if(pts_trk[t][0]->GetN()>0){
                    l_top->AddEntry(pts_trk[t][0],sname.c_str(),"f");
                }else if(pts_trk[t][1]->GetN()>0){
                    l_top->AddEntry(pts_trk[t][1],sname.c_str(),"f");
                }else if(pts_trk[t][2]->GetN()>0){
                    l_top->AddEntry(pts_trk[t][2],sname.c_str(),"f");
                }
            }
            l_top->SetLineWidth(0);
            l_top->SetLineColor(kWhite);
            l_top->Draw("same");

            //Clusters

            /* if(m_is_data==false){ 
               for(auto &c: vec_clusters){
            //auto ssscor = this->ScoreCluster(c.getPlane(),c.getID(), c.getHits() ,vertex_wire[c.getPlane()], vertex_time[c.getPlane()], showers[0]);
            //c.setSSScore(ssscor);

            int thisid = m_sim_shower_trackID[0];

            for(auto &h: c.getHits()){


            }
            }
            }
            */

//CHECK
//            can->cd(8);
//            for(int i=0; i<3; i++){
//                TLegend * l_bot = new TLegend(0.1+i*0.25,0.1,0.1+i*0.25+0.25,0.89);
//
//                TLegend * l_bot2 = new TLegend(0.1+i*0.25,0.1,0.1+i*0.25+0.25,0.89);
//
//                for(int c=0; c< num_clusters[i]+1; c++){
//                    if(c==0)continue;
//
//                    int num_hits_in_cluster = v_newClusterToHitsMap[i][c].size();
//                    auto hitz = v_newClusterToHitsMap[i][c];
//                    auto ssscorz = this->ScoreCluster(i,c, hitz ,vertex_wire[i], vertex_time[i], showers[0]);
//                    int is_in_shower = this->CompareToShowers(i,c, hitz ,vertex_wire[i], vertex_time[i], showers, showerToPFParticleMap, pfParticleToHitsMap,eps);
//
//                    double mean_summed_ADC = 0.0;
//                    for(auto &h:hitz){
//                        mean_summed_ADC +=h->SummedADC();
//                    }
//                    mean_summed_ADC = mean_summed_ADC/(double)num_hits_in_cluster;
//
//
//
//                    // std::string sname = makeSplitlineString({"Cluster: ","Hits: ","PCA: ","Theta: "},{c,num_hits_in_cluster});
//
//                    std::string sname = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Hits: "+std::to_string(num_hits_in_cluster)+"}{#splitline{PCA "+std::to_string(ssscorz.pca_0)+"}{#splitline{Theta:" +std::to_string(ssscorz.pca_theta)+"}{#splitline{Wires: "+std::to_string(ssscorz.n_wires)+ "}{#splitline{Ticks: "+std::to_string(ssscorz.n_ticks)+"}{#splitline{ReMerged: "+std::to_string(is_in_shower)+"}{}}}}}}}";
//                    l_bot->AddEntry(g_clusters[i][c],sname.c_str(),"f");
//
//                    //Here we will only plot those that pass in bottom:
//                    //We are also going to put a hard threshold of 70cm? 
//                    //
//                    if(ssscorz.pass && is_in_shower ==-1 ){
//                        can->cd(i+9);
//                        if(g_clusters[i][c]->GetN()>0){
//                            TGraph * tmp = (TGraph*)g_clusters[i][c]->Clone(("tmp_"+std::to_string(i)+std::to_string(c)).c_str());
//
//                            int Npts =m_SEAviewMaxPtsLinFit;
//                            TGraph * core  = (TGraph*)this->GetNearestNpts(i,c,hitz,vertex_wire[i],vertex_time[i],Npts);
//
//                            core->Draw("p same");
//                            tmp->Draw("p same");
//
//                            double fmax = -999;
//                            double fmin = 99999;
//                            for(int b=0; b<core->GetN(); b++){
//                                double ttx=0;
//                                double tty=0;
//                                core->GetPoint(b,ttx,tty);
//                                fmax = std::max(fmax, ttx);
//                                fmin = std::min(fmin,ttx);
//                            }
//
////                            std::cout<<"Just Before Core Fit of "<<tmp->GetN()<<" pts between "<<chan_min[i]<<" "<<chan_max[i]<<" or "<<fmin<<" "<<fmax<<std::endl;
//
//                            double con;
//                            double slope;
//                            if(fmin==fmax){
//                                slope = 0;
//                                con = fmin;
//                            }else{
//                                std::cout<<sname<<std::endl;
//                                core->Fit("pol1","Q","same",fmin,fmax);
//                                core->GetFunction("pol1")->SetLineWidth(1); 
//                                core->GetFunction("pol1")->SetLineStyle(3); 
//                                core->GetFunction("pol1")->SetLineColor(g_clusters[i][c]->GetMarkerColor()); 
//                                con = core->GetFunction("pol1")->GetParameter(0);
//                                slope = core->GetFunction("pol1")->GetParameter(1);
//                            }
//                            //lets map (wire,tick) to a rudamentary (cm,cm);
//                            //double slope2 = slope*25*0.3;
//                            //double con2 = con*25;
//
//                            double impact_parameter = 1e10;// fabs(slope*vertex_wire[i] +vertex_time[i]+con)/sqrt(slope*slope+1.0*1.0);
//
//                            //rudimentary!
//                            for(double k=chan_min[i]; k< chan_max[i];k++){
//                                double y = slope*k+con;
//                                double dist = sqrt(pow(k*0.3-vertex_wire[i]*0.3,2)+pow(y/25.0-vertex_time[i]/25.0,2));
//                                impact_parameter = std::min(impact_parameter,dist);
//                            }
//
//                            //Lets assume its potining back to the "vertex" and calculate a kinda_angle w.r.t to the shower
//                            //vertex_wire[i] vertex_tick[i] (already calcuated)
//                            //cluster closest point )ssscorz.close_wire and close_tick
//                            //recob::Shower start point, convered to wire tick.
//
//                            double shr_wire = (double)calcWire(m_reco_shower_starty[0], m_reco_shower_startz[0], i, fTPC, fCryostat, *geom);
//                            double shr_time = calcTime(m_reco_shower_startx[0], i, fTPC,fCryostat, *theDetector);
//
//                            std::vector<double> vec_c = {(double)(vertex_wire[i]-ssscorz.close_wire), (double)(vertex_time[i]-ssscorz.close_tick)};
//                            std::vector<double> vec_s = {(double)vertex_wire[i]-shr_wire, (double)vertex_time[i]-shr_time};
//                            double l_c = sqrt(pow(0.3*vec_c[0],2)+pow(vec_c[1]/25.0,2));
//                            double l_s = sqrt(pow(0.3*vec_s[0],2)+pow(vec_s[1]/25.0,2));
//                            double kinda_angle = acos((0.3*vec_s[0]*0.3*vec_c[0]+vec_c[1]*vec_s[1]/(25.0*25.0) )/(l_c*l_s));
//                            //std::cout<<"KINDA "<<kinda_angle<<" "<<l_c<<" "<<l_s<<std::endl;
//
//                            m_sss_num_candidates++;
//
//                            m_sss_candidate_num_hits.push_back(num_hits_in_cluster);
//                            m_sss_candidate_num_wires.push_back((int)ssscorz.n_wires);
//                            m_sss_candidate_num_ticks.push_back((int)ssscorz.n_ticks);
//                            m_sss_candidate_plane.push_back((int)i);
//                            m_sss_candidate_PCA.push_back(ssscorz.pca_0);
//                            m_sss_candidate_impact_parameter.push_back(impact_parameter);
//                            m_sss_candidate_fit_slope.push_back(slope);
//                            m_sss_candidate_fit_constant.push_back(con);
//                            m_sss_candidate_mean_tick.push_back(ssscorz.mean_tick);
//                            m_sss_candidate_max_tick.push_back(ssscorz.max_tick);
//                            m_sss_candidate_min_tick.push_back(ssscorz.min_tick);
//                            m_sss_candidate_min_wire.push_back(ssscorz.min_wire);
//                            m_sss_candidate_max_wire.push_back(ssscorz.max_wire);
//                            m_sss_candidate_mean_wire.push_back(ssscorz.mean_wire);
//                            m_sss_candidate_min_dist.push_back(ssscorz.min_dist);
//                            m_sss_candidate_mean_ADC.push_back(mean_summed_ADC);
//
//                            m_sss_candidate_energy.push_back( this->CalcEShowerPlane(hitz,(int)i));
//                            m_sss_candidate_angle_to_shower.push_back(kinda_angle);
//
//
//                            if(m_is_data){
//                                m_sss_candidate_matched.push_back(-1);
//                                m_sss_candidate_pdg.push_back(-1);
//                                m_sss_candidate_parent_pdg.push_back(-1);
//                                m_sss_candidate_trackid.push_back(-1);
//                                m_sss_candidate_overlay_fraction.push_back(-1);
//				m_sss_candidate_matched_energy_fraction_best_plane.push_back(-1);
//                            }else{
//                                auto ssmatched = this->SecondShowerMatching(hitz, mcparticles_per_hit, mcParticleVector, pfParticleIdMap,  MCParticleToTrackIdMap);
//                                m_sss_candidate_matched.push_back(ssmatched[0]);
//                                m_sss_candidate_pdg.push_back(ssmatched[1]);
//                                m_sss_candidate_parent_pdg.push_back(ssmatched[2]);
//                                m_sss_candidate_trackid.push_back(ssmatched[3]);
//                                m_sss_candidate_overlay_fraction.push_back(ssmatched[4]);
//				m_sss_candidate_matched_energy_fraction_best_plane.push_back(ssmatched[5]);
//                            }
//
//
//                            std::string sname2 = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Impact: "+std::to_string(impact_parameter)+"}{MinDist: "+std::to_string(ssscorz.min_dist)+"}}";
//                            l_bot2->AddEntry(tmp,sname2.c_str(),"f");
//                        }
//                    }          
//
//                }
//
//                //Some time matching
//
//
//                //Closest neightor
//                for(int l=0; l< m_sss_num_candidates; l++){
//                    int this_p = m_sss_candidate_plane[l];
//                    double close = 1e10;                        
//                    for(int m=0; m< m_sss_num_candidates;m++){
//                        if(this_p == m_sss_candidate_plane[m]) continue;
//                        double dup = fabs(m_sss_candidate_mean_tick[l] - m_sss_candidate_mean_tick[m]);
//                        close = std::min(dup,close);
//                    }
//                    m_sss_candidate_closest_neighbour.push_back(close);
//                }
//
//                for(int l=0; l< m_sss_num_candidates; l++){
//
//                    std::vector<double> thisvars = { (double)m_sss_candidate_num_hits[l], (double)m_sss_candidate_num_wires[l], (double)m_sss_candidate_num_ticks[l], (double)m_sss_candidate_PCA[l], log10((double)m_sss_candidate_impact_parameter[l]), log10((double)m_sss_candidate_min_dist[l]), (double)m_sss_candidate_impact_parameter[l]/(double)m_sss_candidate_min_dist[l], (double)m_sss_candidate_energy[l]*0.001, cos((double)m_sss_candidate_angle_to_shower[l]), (double)m_sss_candidate_fit_slope[l], (double)m_sss_candidate_fit_constant[l], (double)m_sss_candidate_plane[l],m_reco_shower_energy_max[0]*0.001,  2*0.001*0.001*m_sss_candidate_energy[l]*m_reco_shower_energy_max[0]*(1.0-cos(m_sss_candidate_angle_to_shower[l])) , log10(2*0.001*0.001*m_sss_candidate_energy[l]*m_reco_shower_energy_max[0]*(1.0-cos(m_sss_candidate_angle_to_shower[l]))),m_sss_candidate_energy[l]/m_reco_shower_energy_max[0], (double)m_sss_candidate_closest_neighbour[l] };
//
//
//                    //double score = sssVetov1->GetMvaValue(thisvars);
//                    double score = -1;
//                    m_sss_candidate_veto_score.push_back(score);
//
//                }
//
//
//                can->cd(8);
//                l_bot->SetLineWidth(0);
//                l_bot->SetLineColor(kWhite);
//                l_bot->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
//                l_bot->Draw("same");
//
//                can->cd(12);
//                l_bot2->SetLineWidth(0);
//                l_bot2->SetLineColor(kWhite);
//                l_bot2->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
//                l_bot2->Draw("same");
//
//
//
//
//
//
//
//
//            }
            //********** Some Error Checking ********************//

            /*for(int i=0; i<3; i++){

              std::cout<<"Plane "<<i<<" Vertex pts "<<g_vertex[i]->GetN()<<std::endl;
              for(size_t s=0; s< pts_shr.size(); s++){
              std::cout<<"Plane "<<i<<" Shower "<<s<<" pts "<<pts_shr[s][i]->GetN()<<std::endl;
              }
              for(size_t t=0;t < pts_trk.size(); t++){
              std::cout<<"Plane "<<i<<" Track "<<t<<" pts "<<pts_trk[t][i]->GetN()<<std::endl;
              }
              for(int c=0; c< num_clusters[i]+1; c++){
              std::cout<<"Plane "<<i<<" Cluster "<<c<<" pts "<<g_clusters[i][c]->GetN()<<std::endl;
              }
              }*/





            if(false){
                std::cout<<"SinglePhoton::SSS\t||\tDone Plotting clusters"<<std::endl;
                can->Update();
                //can->Write();
                can->SaveAs((print_name+".pdf").c_str(),"pdf");
                //f->Close();
                std::cout<<"SinglePhoton::SSS\t||\tPRINTING"<<std::endl;
            }
            //bool_make_sss_plots=false;
            delete can;






            }

        }else{

            std::cout<<"SinglePhoton::SSS\t||\tNo Neutrino Slice Found/No Showers & tracks"<<std::endl;

        }


        //For purposes of testing SEAview compatability:
        //

        std::cout<<"SSSOLD "<<m_run_number<<" "<<m_subrun_number<<" "<<m_event_number<<" NumCandidates: "<<m_sss_num_candidates<<std::endl;
         for(int i=0; i< m_sss_num_candidates; i++){
            std::cout<<i<<" Num Hits "<<m_sss_candidate_num_hits[i]<<std::endl;
            std::cout<<i<<" Num Wires "<<m_sss_candidate_num_wires[i]<<std::endl;
            std::cout<<i<<" Num Ticks "<<m_sss_candidate_num_ticks[i]<<std::endl;
            std::cout<<i<<" Plane "<<m_sss_candidate_plane[i]<<std::endl;
            std::cout<<i<<" PCA "<<m_sss_candidate_PCA[i]<<std::endl;
            std::cout<<i<<" Impact "<<m_sss_candidate_impact_parameter[i]<<std::endl;
            std::cout<<i<<" Fit Slope "<<m_sss_candidate_fit_slope[i]<<std::endl;
            std::cout<<i<<" Fit Constant "<<m_sss_candidate_fit_constant[i]<<std::endl;
            std::cout<<i<<" Mean Tick "<<m_sss_candidate_mean_tick[i]<<std::endl;
            std::cout<<i<<" Max Tick "<<m_sss_candidate_max_tick[i]<<std::endl;
            std::cout<<i<<" Min Tick "<<m_sss_candidate_min_tick[i]<<std::endl;
            std::cout<<i<<" Mean Wire "<<m_sss_candidate_mean_wire[i]<<std::endl;
            std::cout<<i<<" Max Wire "<<m_sss_candidate_max_wire[i]<<std::endl;
            std::cout<<i<<" Min Wire "<<m_sss_candidate_min_wire[i]<<std::endl;
            std::cout<<i<<" Energy "<<m_sss_candidate_energy[i]<<std::endl;
            std::cout<<i<<" AngletoShower "<<m_sss_candidate_angle_to_shower[i]<<std::endl;
        }

        return;
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

        std::string sss3dlabel = "allShr";//"pandoraAllOutcomesShower"
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
