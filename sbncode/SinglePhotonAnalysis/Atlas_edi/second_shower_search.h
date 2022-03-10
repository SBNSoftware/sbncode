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
        m_sss_num_associated_hits=0;

        m_sss_num_candidates = 0;

        m_sss_candidate_num_hits.clear();
        m_sss_candidate_num_wires.clear();
        m_sss_candidate_num_ticks.clear();
        m_sss_candidate_plane.clear();
        m_sss_candidate_PCA.clear();
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
        m_sss_candidate_energy.clear();
        m_sss_candidate_angle_to_shower.clear();
        m_sss_candidate_closest_neighbour.clear();
        m_sss_candidate_matched.clear();
        m_sss_candidate_pdg.clear();
        m_sss_candidate_parent_pdg.clear();
        m_sss_candidate_trackid.clear();
        m_sss_candidate_overlay_fraction.clear();
    }

    void SinglePhoton::ResizeSecondShowers(size_t size){

    }


    void SinglePhoton::CreateSecondShowerBranches(){
        vertex_tree->Branch("sss_num_unassociated_hits",&m_sss_num_unassociated_hits,"sss_num_unassociated_hits/I");
        vertex_tree->Branch("sss_num_associated_hits",&m_sss_num_associated_hits,"sss_num_associated_hits/I");

        vertex_tree->Branch("sss_num_candidates",&m_sss_num_candidates,"sss_num_candidates/I");
        vertex_tree->Branch("sss_candidate_veto_score",&m_sss_candidate_veto_score);
        vertex_tree->Branch("sss_candidate_num_hits",&m_sss_candidate_num_hits);
        vertex_tree->Branch("sss_candidate_num_wires",&m_sss_candidate_num_wires);
        vertex_tree->Branch("sss_candidate_num_ticks",&m_sss_candidate_num_ticks);
        vertex_tree->Branch("sss_candidate_plane",&m_sss_candidate_plane);
        vertex_tree->Branch("sss_candidate_PCA",&m_sss_candidate_PCA);
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
        vertex_tree->Branch("sss_candidate_energy",&m_sss_candidate_energy);
        vertex_tree->Branch("sss_candidate_angle_to_shower",&m_sss_candidate_angle_to_shower);
        vertex_tree->Branch("sss_candidate_closest_neighbour",&m_sss_candidate_closest_neighbour);

        vertex_tree->Branch("sss_candidate_matched",&m_sss_candidate_matched);
        vertex_tree->Branch("sss_candidate_pdg",&m_sss_candidate_pdg);
        vertex_tree->Branch("sss_candidate_parent_pdg",&m_sss_candidate_parent_pdg);
        vertex_tree->Branch("sss_candidate_trackid",&m_sss_candidate_trackid);
        vertex_tree->Branch("sss_candidate_overlay_fraction",&m_sss_candidate_overlay_fraction);


    }

    void SinglePhoton::SecondShowerSearch(
            const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
            const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
            const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
            const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
            std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap){



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
		std::cout<<"CHECK1! in second_shower_search.h"<<endl;
            auto slicehits = sliceIDToHitsMap.at(sliceid);
		std::cout<<"CHECK2!"<<endl;
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


            std::vector<TGraph*> g_vertex(3);//the star mark
            for(int i=0; i<3; i++){
                TPad * pader = (TPad*)can->cd(i+1);

                if(i==0 || i ==4 || i == 8) pader->SetLeftMargin(0.1);

                std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *geom)};
                std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, *theDetector)};

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
            int min_pts = 10;
            double eps = 5.0;
            std::vector<int> num_clusters(3,0);

            std::vector<std::vector<TGraph*>> g_clusters(3);
            std::vector<std::vector<int>> cluster_labels(3);


            std::vector<cluster> vec_clusters;

            for(int i=0; i<3; i++){

                std::cout<<"SinglePhoton::SSS\t||\tStarting to run DBSCAN for plane: "<<i<<" has "<<pts_to_recluster[i].size()<<" pts to do using eps: "<<eps<<" and min_pts: "<<min_pts<<std::endl; 
                DBSCAN ReCluster(eps,min_pts);
                cluster_labels[i] =  ReCluster.Scan2D(pts_to_recluster[i]);

                for(auto &c: cluster_labels[i]){
                    num_clusters[i] = std::max(c,num_clusters[i]);
                }
                std::cout << "SinglePhoton::SSS\t||\tOn this plane, DBSCAN found: "<<num_clusters[i]<<" clusters"<<std::endl;
            }

            //Now fill the clusters

            for(int i=0; i<3; i++){

                for(int c=0; c<num_clusters[i]; c++){

                    std::vector<std::vector<double>> pts;
                    std::vector<art::Ptr<recob::Hit>> hitz;
                    for(size_t p=0; p< pts_to_recluster[i].size(); p++){
                        if(cluster_labels[i][p] == 0) continue;//noise 
                        if(cluster_labels[i][p] == c){

                            pts.push_back(pts_to_recluster[i][p]);
                            hitz.push_back(mapPointIndexToHit[i].at(p));
                        }

                    }
                    vec_clusters.emplace_back(c,i,pts,hitz);
                }
            }



            //Create final maps;
            for(size_t i=0; i<3; i++){
                for(size_t p=0; p<pts_to_recluster[i].size(); p++){

                    art::Ptr<recob::Hit> h = mapPointIndexToHit[i].at(p);// Get the hit
                    size_t cluster_label = cluster_labels[i][p];//get the cluster index, 0 = noise                                                

                    //std::cout<<i<<" "<<p<<" "<<cluster_label<<" "<<v_newClusterToHitsMap[i].count(cluster_label)<<std::endl;

                    if(v_newClusterToHitsMap[i].count(cluster_label)<1){
                        std::vector<art::Ptr<recob::Hit>> t_hs= {h};
                        v_newClusterToHitsMap[i][cluster_label] = t_hs;
                    }else{
                        v_newClusterToHitsMap[i].at(cluster_label).push_back(h);//add it to list
                    }
                }
            }


            //for(size_t i=0; i<3; i++){
            //  for(int c=0; c<num_clusters[i]+1; c++){
            //auto v_newClusterToHitsMap[i][c];
            // }
            // }


            int max_n_clusters = std::max(num_clusters[0]+1, std::max(num_clusters[1]+1,num_clusters[2]+1))+2;
            std::vector<int> cluster_colors(max_n_clusters,0);
            std::vector<int> base_col = {632,416, 600, 400, 616,  432,  800, 820,  840, 860,  880, 900};

            for(int j=0; j< max_n_clusters; j++){
                int b = (int)rangen->Uniform(0,11);
                int mod = (int)rangen->Uniform(-10,+3);

                cluster_colors[j] = base_col[b]+mod;
            }
            int c_offset = 0;

            //Step next, loop over and make plots again
            for(int i=0; i<3; i++){
                std::vector<std::vector<double>> vec_time(num_clusters[i]+1);
                std::vector<std::vector<double>> vec_wire(num_clusters[i]+1);
                std::vector<TGraph*> tmp_g_clusters(num_clusters[i]+1);

                if(cluster_labels[i].size() != pts_to_recluster[i].size()){
                    std::cout<<"SinglePhoton::SSS\t||\tERROR!! someting amiss cluster labels of size "<<cluster_labels[i].size()<<" and pts  in this plane "<<pts_to_recluster[i].size()<<std::endl;
                }  

                for(size_t k=0; k< pts_to_recluster[i].size(); k++){
                    //std::cout<<vec_wire.size()<<" "<<cluster_labels[i][k]<<std::endl;
                    vec_wire[cluster_labels[i][k]].push_back(pts_to_recluster[i][k][0]); 
                    vec_time[cluster_labels[i][k]].push_back(pts_to_recluster[i][k][1]); 

                }

                for(int c=0; c< num_clusters[i]+1; c++){
                    int tcol = kBlack;
                    if(c>0) tcol = cluster_colors[c+c_offset];
                    tmp_g_clusters[c] = new TGraph(vec_wire[c].size(),&(vec_wire[c])[0],&(vec_time[c])[0] );
                    can->cd(i+5);
                    if(
                            tmp_g_clusters[c]->GetN()>0){
                        tmp_g_clusters[c]->Draw("p same");
                        tmp_g_clusters[c]->SetMarkerColor(tcol);
                        tmp_g_clusters[c]->SetFillColor(tcol);
                        tmp_g_clusters[c]->SetMarkerStyle(20);
                        tmp_g_clusters[c]->SetMarkerSize(plot_point_size);
                    }
                }
                g_clusters[i] = tmp_g_clusters;
                c_offset += num_clusters[i];
            }

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


            can->cd(8);
            for(int i=0; i<3; i++){
                TLegend * l_bot = new TLegend(0.1+i*0.25,0.1,0.1+i*0.25+0.25,0.89);

                TLegend * l_bot2 = new TLegend(0.1+i*0.25,0.1,0.1+i*0.25+0.25,0.89);

                for(int c=0; c< num_clusters[i]+1; c++){
                    if(c==0)continue;

                    int num_hits_in_cluster = v_newClusterToHitsMap[i][c].size();
                    auto hitz = v_newClusterToHitsMap[i][c];
                    auto ssscorz = this->ScoreCluster(i,c, hitz ,vertex_wire[i], vertex_time[i], showers[0]);
                    int is_in_shower = this->CompareToShowers(i,c, hitz ,vertex_wire[i], vertex_time[i], showers, showerToPFParticleMap, pfParticleToHitsMap,eps);

                    // std::string sname = makeSplitlineString({"Cluster: ","Hits: ","PCA: ","Theta: "},{c,num_hits_in_cluster});

                    std::string sname = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Hits: "+std::to_string(num_hits_in_cluster)+"}{#splitline{PCA "+std::to_string(ssscorz.pca_0)+"}{#splitline{Theta:" +std::to_string(ssscorz.pca_theta)+"}{#splitline{Wires: "+std::to_string(ssscorz.n_wires)+ "}{#splitline{Ticks: "+std::to_string(ssscorz.n_ticks)+"}{#splitline{ReMerged: "+std::to_string(is_in_shower)+"}{}}}}}}}";
                    l_bot->AddEntry(g_clusters[i][c],sname.c_str(),"f");

                    //Here we will only plot those that pass in bottom:
                    //We are also going to put a hard threshold of 70cm? 
                    //
                    if(ssscorz.pass && is_in_shower ==-1 ){
                        can->cd(i+9);
                        if(g_clusters[i][c]->GetN()>0){
                            TGraph * tmp = (TGraph*)g_clusters[i][c]->Clone(("tmp_"+std::to_string(i)+std::to_string(c)).c_str());

                            int Npts = 20;
                            TGraph * core  = (TGraph*)this->GetNearestNpts(i,c,hitz,vertex_wire[i],vertex_time[i],Npts);

                            core->Draw("p same");
                            tmp->Draw("p same");

                            double fmax = -999;
                            double fmin = 99999;
                            for(int b=0; b<core->GetN(); b++){
                                double ttx=0;
                                double tty=0;
                                core->GetPoint(b,ttx,tty);
                                fmax = std::max(fmax, ttx);
                                fmin = std::min(fmin,ttx);
                            }

                            std::cout<<"Just Before Core Fit of "<<tmp->GetN()<<" pts between "<<chan_min[i]<<" "<<chan_max[i]<<" or "<<fmin<<" "<<fmax<<std::endl;

                            double con;
                            double slope;
                            if(fmin==fmax){
                                slope = 0;
                                con = fmin;
                            }else{
                                std::cout<<sname<<std::endl;
                                core->Fit("pol1","Q","same",fmin,fmax);
                                core->GetFunction("pol1")->SetLineWidth(1); 
                                core->GetFunction("pol1")->SetLineStyle(3); 
                                core->GetFunction("pol1")->SetLineColor(g_clusters[i][c]->GetMarkerColor()); 
                                con = core->GetFunction("pol1")->GetParameter(0);
                                slope = core->GetFunction("pol1")->GetParameter(1);

                            }
                            //lets map (wire,tick) to a rudamentary (cm,cm);
                            //double slope2 = slope*25*0.3;
                            //double con2 = con*25;

                            double impact_parameter = 1e10;// fabs(slope*vertex_wire[i] +vertex_time[i]+con)/sqrt(slope*slope+1.0*1.0);

                            //rudimentary!
                            for(double k=chan_min[i]; k< chan_max[i];k++){
                                double y = slope*k+con;
                                double dist = sqrt(pow(k*0.3-vertex_wire[i]*0.3,2)+pow(y/25.0-vertex_time[i]/25.0,2));
                                impact_parameter = std::min(impact_parameter,dist);
                            }


                            //Lets assume its potining back to the "vertex" and calculate a kinda_angle w.r.t to the shower
                            //vertex_wire[i] vertex_tick[i] (already calcuated)
                            //cluster closest point )ssscorz.close_wire and close_tick
                            //recob::Shower start point, convered to wire tick.

                            double shr_wire = (double)calcWire(m_reco_shower_starty[0], m_reco_shower_startz[0], i, fTPC, fCryostat, *geom);
                            double shr_time = calcTime(m_reco_shower_startx[0], i, fTPC,fCryostat, *theDetector);

                            std::vector<double> vec_c = {(double)(vertex_wire[i]-ssscorz.close_wire), (double)(vertex_time[i]-ssscorz.close_tick)};
                            std::vector<double> vec_s = {(double)vertex_wire[i]-shr_wire, (double)vertex_time[i]-shr_time};
                            double l_c = sqrt(pow(0.3*vec_c[0],2)+pow(vec_c[1]/25.0,2));
                            double l_s = sqrt(pow(0.3*vec_s[0],2)+pow(vec_s[1]/25.0,2));
                            double kinda_angle = acos((0.3*vec_s[0]*0.3*vec_c[0]+vec_c[1]*vec_s[1]/(25.0*25.0) )/(l_c*l_s));
                            //std::cout<<"KINDA "<<kinda_angle<<" "<<l_c<<" "<<l_s<<std::endl;

                            m_sss_num_candidates++;

                            m_sss_candidate_num_hits.push_back(num_hits_in_cluster);
                            m_sss_candidate_num_wires.push_back((int)ssscorz.n_wires);
                            m_sss_candidate_num_ticks.push_back((int)ssscorz.n_ticks);
                            m_sss_candidate_plane.push_back((int)i);
                            m_sss_candidate_PCA.push_back(ssscorz.pca_0);
                            m_sss_candidate_impact_parameter.push_back(impact_parameter);
                            m_sss_candidate_fit_slope.push_back(slope);
                            m_sss_candidate_fit_constant.push_back(con);
                            m_sss_candidate_mean_tick.push_back(ssscorz.mean_tick);
                            m_sss_candidate_max_tick.push_back(ssscorz.max_tick);
                            m_sss_candidate_min_tick.push_back(ssscorz.min_tick);
                            m_sss_candidate_min_wire.push_back(ssscorz.min_wire);
                            m_sss_candidate_max_wire.push_back(ssscorz.max_wire);
                            m_sss_candidate_mean_wire.push_back(ssscorz.mean_wire);
                            m_sss_candidate_min_dist.push_back(ssscorz.min_dist);

                            m_sss_candidate_energy.push_back(  this->CalcEShowerPlane(hitz,(int)i));
                            m_sss_candidate_angle_to_shower.push_back(kinda_angle);


                            if(m_is_data){
                                m_sss_candidate_matched.push_back(-1);
                                m_sss_candidate_pdg.push_back(-1);
                                m_sss_candidate_parent_pdg.push_back(-1);
                                m_sss_candidate_trackid.push_back(-1);
                                m_sss_candidate_overlay_fraction.push_back(-1);

                            }else{
                                auto ssmatched = this->SecondShowerMatching( hitz, mcparticles_per_hit, mcParticleVector, pfParticleIdMap,  MCParticleToTrackIdMap);
                                m_sss_candidate_matched.push_back(ssmatched[0]);
                                m_sss_candidate_pdg.push_back(ssmatched[1]);
                                m_sss_candidate_parent_pdg.push_back(ssmatched[2]);
                                m_sss_candidate_trackid.push_back(ssmatched[3]);
                                m_sss_candidate_overlay_fraction.push_back(ssmatched[4]);


                            }


                            std::string sname2 = "#splitline{Cluster "+std::to_string(c)+"}{#splitline{Impact: "+std::to_string(impact_parameter)+"}{MinDist: "+std::to_string(ssscorz.min_dist)+"}}";
                            l_bot2->AddEntry(tmp,sname2.c_str(),"f");
                        }
                    }          

                }

                //Some time matching

                //Closest neightor

                for(int l=0; l< m_sss_num_candidates; l++){
                    int this_p = m_sss_candidate_plane[l];
                    double close = 1e10;                        
                    for(int m=0; m< m_sss_num_candidates;m++){
                        if(this_p == m_sss_candidate_plane[m]) continue;

                        double dup = fabs(m_sss_candidate_mean_tick[l] - m_sss_candidate_mean_tick[m]);

                        close = std::min(dup,close);

                    }

                    m_sss_candidate_closest_neighbour.push_back(close);

                }



                for(int l=0; l< m_sss_num_candidates; l++){

                    std::vector<double> thisvars = { (double)m_sss_candidate_num_hits[l], (double)m_sss_candidate_num_wires[l], (double)m_sss_candidate_num_ticks[l], (double)m_sss_candidate_PCA[l], log10((double)m_sss_candidate_impact_parameter[l]), log10((double)m_sss_candidate_min_dist[l]), (double)m_sss_candidate_impact_parameter[l]/(double)m_sss_candidate_min_dist[l], (double)m_sss_candidate_energy[l]*0.001, cos((double)m_sss_candidate_angle_to_shower[l]), (double)m_sss_candidate_fit_slope[l], (double)m_sss_candidate_fit_constant[l], (double)m_sss_candidate_plane[l],m_reco_shower_energy_max[0]*0.001,  2*0.001*0.001*m_sss_candidate_energy[l]*m_reco_shower_energy_max[0]*(1.0-cos(m_sss_candidate_angle_to_shower[l])) , log10(2*0.001*0.001*m_sss_candidate_energy[l]*m_reco_shower_energy_max[0]*(1.0-cos(m_sss_candidate_angle_to_shower[l]))),m_sss_candidate_energy[l]/m_reco_shower_energy_max[0], (double)m_sss_candidate_closest_neighbour[l] };
                    

                    double score = sssVetov1->GetMvaValue(thisvars);
                    m_sss_candidate_veto_score.push_back(score);

                }


                can->cd(8);
                l_bot->SetLineWidth(0);
                l_bot->SetLineColor(kWhite);
                l_bot->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
                l_bot->Draw("same");

                can->cd(12);
                l_bot2->SetLineWidth(0);
                l_bot2->SetLineColor(kWhite);
                l_bot2->SetHeader(("Plane "+std::to_string(i)).c_str(),"C");
                l_bot2->Draw("same");








            }
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


        return;
    }



    TGraph* SinglePhoton::GetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts){

        std::vector<double>t_wire;
        std::vector<double>t_tick;
        // std::vector<double>t_dist;

        std::vector<double>all_wire;
        std::vector<double>all_tick;
        std::vector<double>all_dist;


        for(size_t h = 0; h< hitz.size(); h++){
            auto hit = hitz[h];
            double h_wire = (double)hit->WireID().Wire;
            double h_tick = (double)hit->PeakTime();

            double dd =sqrt(pow(h_wire*0.3-vertex_wire*0.3,2)+pow(h_tick/25.0- vertex_tick/25.0,2));
            all_wire.push_back(h_wire);   
            all_tick.push_back(h_tick);   
            all_dist.push_back(dd);
        }

        std::vector<size_t> sorted_in = sort_indexes(all_dist);
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



    std::vector<double>SinglePhoton::SecondShowerMatching(std::vector<art::Ptr<recob::Hit>>& hitz,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
            std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap){


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


        }
        if(found_a_match){
            std::cout<<"!"<<std::endl;
        }
        double fraction_num_hits_overlay = (double)n_not_associated_hits/(double)hitz.size();

        //            if(reco_verbose)std::cout << "SinglePhoton::recoMC()\t||\t On Object "<<i<<". The number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       
        //          if(reco_verbose) std::cout<<"SinglePhoton::recoMC()\t||\t the fraction of hits from overlay is is "<<fraction_num_hits_overlay<<" ("<<n_not_associated_hits<<"/"<<obj_hits_ptrs.size()<<")"<<std::endl;


        if(n_associated_mcparticle_hits == 0){
            //This will only occur if the whole recob::PFParticle is PURELY associated with an overlay object
            found_a_match =false;
            //Here we will fill every sim_shower_XXX variable with -999 or something like that 
            return {0,0,0,0,0};
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

        ans = {1,(double)match->PdgCode(), (double)par_pdg, (double)match->TrackId(),fraction_num_hits_overlay};

        return ans;
    }//end sss matching;

}
