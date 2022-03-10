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

#include "TH1.h"

// override function of sorts for max_element function comparison
bool  map_max_fn(const std::pair<art::Ptr<recob::Hit>,double> p1, const std::pair<art::Ptr<recob::Hit>,  double> p2){
	return (p1.second < p2.second);
}

// override function of sorts for min_element function comparison
bool  map_min_fn(const std::pair<art::Ptr<recob::Hit>,double> p1, const std::pair<art::Ptr<recob::Hit>,  double> p2){
	return (p1.second > p2.second);
}

namespace single_photon{

void SinglePhoton::ClearIsolation(){
    m_isolation_min_dist_trk_shr.clear();
    m_isolation_min_dist_trk_unassoc.clear();

    m_isolation_num_shr_hits_win_1cm_trk.clear();
    m_isolation_num_shr_hits_win_2cm_trk.clear();
    m_isolation_num_shr_hits_win_5cm_trk.clear();
    m_isolation_num_shr_hits_win_10cm_trk.clear();

    m_isolation_num_unassoc_hits_win_1cm_trk.clear();
    m_isolation_num_unassoc_hits_win_2cm_trk.clear();
    m_isolation_num_unassoc_hits_win_5cm_trk.clear();
    m_isolation_num_unassoc_hits_win_10cm_trk.clear();

    m_isolation_nearest_shr_hit_to_trk_wire.clear();
    m_isolation_nearest_shr_hit_to_trk_time.clear();

    m_isolation_nearest_unassoc_hit_to_trk_wire.clear();
    m_isolation_nearest_unassoc_hit_to_trk_time.clear();
}

void SinglePhoton::CreateIsolationBranches(){
    vertex_tree->Branch("isolation_min_dist_trk_shr", &m_isolation_min_dist_trk_shr);
    vertex_tree->Branch("isolation_min_dist_trk_unassoc", &m_isolation_min_dist_trk_unassoc);

    vertex_tree->Branch("isolation_num_shr_hits_win_1cm_trk", &m_isolation_num_shr_hits_win_1cm_trk);
    vertex_tree->Branch("isolation_num_shr_hits_win_2cm_trk", &m_isolation_num_shr_hits_win_2cm_trk);
    vertex_tree->Branch("isolation_num_shr_hits_win_5cm_trk", &m_isolation_num_shr_hits_win_5cm_trk);
    vertex_tree->Branch("isolation_num_shr_hits_win_10cm_trk", &m_isolation_num_shr_hits_win_10cm_trk);
    
    vertex_tree->Branch("isolation_num_unassoc_hits_win_1cm_trk", &m_isolation_num_unassoc_hits_win_1cm_trk);
    vertex_tree->Branch("isolation_num_unassoc_hits_win_2cm_trk", &m_isolation_num_unassoc_hits_win_2cm_trk);
    vertex_tree->Branch("isolation_num_unassoc_hits_win_5cm_trk", &m_isolation_num_unassoc_hits_win_5cm_trk);
    vertex_tree->Branch("isolation_num_unassoc_hits_win_10cm_trk", &m_isolation_num_unassoc_hits_win_10cm_trk);


    vertex_tree->Branch("isolation_nearest_shr_hit_to_trk_wire", &m_isolation_nearest_shr_hit_to_trk_wire);
    vertex_tree->Branch("isolation_nearest_shr_hit_to_trk_time", &m_isolation_nearest_shr_hit_to_trk_time);
    
    vertex_tree->Branch("isolation_nearest_unassoc_hit_to_trk_wire", &m_isolation_nearest_unassoc_hit_to_trk_wire);
    vertex_tree->Branch("isolation_nearest_unassoc_hit_to_trk_time", &m_isolation_nearest_unassoc_hit_to_trk_time);

}


/* Arguments to Function  IsolationStudy  (all are const):
 * 1. vector named tracks 		of art ptr to recob track
 * 2. map named trackToPFPParticleMap 	of .i. art ptr to recob track		.ii. art ptr to recob pfparticle  
 * 3. vector named showers 	      	of art ptr to recob showers
 * 4. map named showerToPFParticleMap 	of .i. art ptr to recob showe		.ii. art ptr to recob pfparticle
 * 5. map named pfParticleToHistMap   	of .i. art ptr to recob prparticle	.ii. vector of art ptr to recob hit
 * 6. map named pfParticleToSliceIDMap	of .i. art ptr to recob pfparticle	.ii. int
 * 7. map named sliceIDToHitsMap	of .i. int and				.ii. art ptr to recob hit
*/
void SinglePhoton::IsolationStudy(
	const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & trackToPFParticleMap,
	const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,
        const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,  
        const std::map<art::Ptr<recob::PFParticle>, int> & pfParticleToSliceIDMap, const std::map<int, std::vector<art::Ptr<recob::Hit>>>& sliceIDToHitsMap)
{

        int total_track_hits =0;
        int total_shower_hits =0;
        int nu_slice_id = -999;

        std::vector< art::Ptr<recob::Hit> > associated_hits;
        std::vector< art::Ptr<recob::Hit> > unassociated_hits;
        std::vector< art::Ptr<recob::Hit> > unassociated_hits_plane0;
        std::vector< art::Ptr<recob::Hit> > unassociated_hits_plane1;
        std::vector< art::Ptr<recob::Hit> > unassociated_hits_plane2;
 
        std::vector< std::map<size_t, std::vector< art::Ptr<recob::Hit> >> > v_newClusterToHitsMap(3);

// BEGIN FOR LOOP TO COUNT TRACK HITS AND PUSH INTO ASSOCIATED HITS
        for(size_t t =0; t< tracks.size(); t++){
	    art::Ptr<recob::Track> track = tracks[t];
            art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];
            
	    int sliceid = pfParticleToSliceIDMap.at(pfp);
            
	    std::vector<art::Ptr<recob::Hit>> slicehits = sliceIDToHitsMap.at(sliceid);
            std::vector<art::Ptr<recob::Hit>> trackhits = pfParticleToHitsMap.at(pfp);

            std::cout << "*SSS: track "<< t <<" is in slice "<< sliceid <<" which has "<< slicehits.size() <<" hits. This track has  "<< trackhits.size() <<" of them. " << std::endl;
            total_track_hits += trackhits.size();

            if(nu_slice_id !=  sliceid && nu_slice_id != -999){
                std::cout<<"*ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
                exit(EXIT_FAILURE);
            }   
            nu_slice_id = sliceid;


            for(auto &h: trackhits){
                associated_hits.push_back(h);
            }

        }
// END FOR LOOPING TRACKS

// BEGIN FOR LOOPING SHOWERS TO COUNT SHOWER HITS AND PUSH INTO ASSOC HITS
        for(size_t s =0; s< showers.size(); s++){
            art::Ptr<recob::Shower> shower = showers[s];
            art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(shower);
            
	    int sliceid = pfParticleToSliceIDMap.at(pfp);
            
	    auto slicehits = sliceIDToHitsMap.at(sliceid);
            auto showerhits = pfParticleToHitsMap.at(pfp);

            std::cout<<"*SSS: shower "<< s <<" is in slice "<< sliceid <<" which has "<< slicehits.size() <<" hits. This shower has  "<< showerhits.size() <<" of them. "<< std::endl;
            total_shower_hits+=showerhits.size();

            if(nu_slice_id !=  sliceid && nu_slice_id!=-999){
                std::cout<<"*ERROR!! In Second Shower Search, the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
                exit(EXIT_FAILURE);
            }
            nu_slice_id = sliceid;

            for(auto &h: showerhits){
                associated_hits.push_back(h);
            } 
        } 
// END FOR LOOP COUNTING SHOWER HITS

// PRINT SUMMARY OF HIT TYPES
        std::cout<<"*SSS: So in total we have "<<total_shower_hits<<" shower hits and "<<total_track_hits<<" track hits"<<" associatedHits total: "<<associated_hits.size()<<std::endl;
        m_sss_num_associated_hits = total_shower_hits + total_track_hits;

// IF VALID SLICE
        if(nu_slice_id >= 0){
            std::cout<<"*SSS: So that leaves "<<sliceIDToHitsMap.at(nu_slice_id).size()-total_shower_hits-total_track_hits<<" hits not included in tracks and showers"<<std::endl;
            m_sss_num_unassociated_hits = sliceIDToHitsMap.at(nu_slice_id).size()-total_shower_hits-total_track_hits;

            if(m_sss_num_unassociated_hits < 0){
                std::cout<<"ERROR!! Number of unassociated hits is negative, i.e: num_associated: "<<m_sss_num_associated_hits<<" and total slice hits: "<<sliceIDToHitsMap.at(nu_slice_id).size()<<std::endl;
                exit(EXIT_FAILURE);
            }



// DETERMINE UNASSOCIATED HITS BY COMPARING ALL SLICE HITS WITH LIST OF ASSOCIATED HITS
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

            std::cout<<" *associated_hits.size() "<<associated_hits.size()<<" unassociated_hits.size() "<<unassociated_hits.size()<<" p0: "<<unassociated_hits_plane0.size()<<" p1:  "<<unassociated_hits_plane1.size()<<" p2: "<<unassociated_hits_plane2.size()<<std::endl;


// IF HAVE 1+G 1P AND WANT TO PLOT
	    if(bool_make_sss_plots && showers.size()==1 && tracks.size()==1){

                std::string print_name = "isolation_"+std::to_string(m_run_number)+"_"+std::to_string(m_subrun_number)+"_"+std::to_string(m_event_number);
                TCanvas *can = new TCanvas(print_name.c_str(), print_name.c_str(),3000,800);
                can->Divide(4, 1, 0.0, 0.1);

                double tick_max = 0;
                double tick_min = 1e10;
                std::vector<double> chan_max(3,0);
                std::vector<double> chan_min(3,1e10);

		// Creation of canvas and histograms to hold hit distance data
		TCanvas *histcan = new TCanvas(("hists_"+print_name).c_str(), "Distances from the track", 600, 400);
		histcan->Divide(3, 2, 0.005, 0.1);
		
		TH1D *s_hist0 = new TH1D("shower hits hist plane0", "Showers Plane 0", 30, 0.0, 30.0);
		TH1D *s_hist1 = new TH1D("shower hist plane1", "Showers Plane 1", 30, 0.0, 30.0);
		TH1D *s_hist2 = new TH1D("shower hist plane2", "Showers Plane 2", 30, 0.0, 30.0);
		std::vector<TH1D *> s_hists = {s_hist0, s_hist1, s_hist2};
		
		TH1D *u_hist0 = new TH1D("unassociated hits hist plane0", "Unassoc Plane 0", 30, 0.0, 30.0);
		TH1D *u_hist1 = new TH1D("unassociated hits hist plane1", "Unassoc Plane 1", 30, 0.0, 30.0);
		TH1D *u_hist2 = new TH1D("unassociated hits hist plane2", "Unassoc Plane 2", 30, 0.0, 30.0);
		std::vector<TH1D *> u_hists = {u_hist0, u_hist1, u_hist2};


std::cout << "Isolation: Acquiring track hit coordinates" << std::endl;
// saving wire and time coordinates
                std::vector<std::vector<TGraph *>> pts_trk( tracks.size(), std::vector<TGraph *>(3)  );
                
                art::Ptr<recob::PFParticle> pfpt = trackToPFParticleMap.at(tracks[0]);
                auto trackhits = pfParticleToHitsMap.at(pfpt);

                std::vector<TGraph*> t_pts(3);
                std::vector<std::vector<double>>  t_vec_t(3);
                std::vector<std::vector<double>>  t_vec_c(3);

                for(auto &th: trackhits){
		    double wire = (double)th->WireID().Wire;
		    t_vec_c[(int)th->View()].push_back(wire);

		    double time = (double)th->PeakTime();
                    t_vec_t[(int)th->View()].push_back(time);

                    tick_max = std::max(tick_max, time);
                    tick_min = std::min(tick_min, time);
                    chan_max[(int)th->View()] = std::max( chan_max[(int)th->View()], wire);
                    chan_min[(int)th->View()] = std::min( chan_min[(int)th->View()], wire);

                }

                t_pts[0] = new TGraph(t_vec_c[0].size(), &(t_vec_c[0])[0], &(t_vec_t[0])[0]);
                t_pts[1] = new TGraph(t_vec_c[1].size(), &(t_vec_c[1])[0], &(t_vec_t[1])[0]);
                t_pts[2] = new TGraph(t_vec_c[2].size(), &(t_vec_c[2])[0], &(t_vec_t[2])[0]);
                pts_trk[0] = t_pts;

std::cout << "Isolation: plane0 track hits = " << t_vec_t[0].size() << std::endl;
std::cout << "Isolation: plane1 track hits = " << t_vec_t[1].size() << std::endl;
std::cout << "Isolation: plane2 track hits = " << t_vec_t[2].size() << std::endl;


std::cout << "Isolation: Acquiring  shower hit coordinates and comparing with track hits " << std::endl;
		
		//First grab all shower clusters
		 std::vector<std::vector<TGraph *>> pts_shr( showers.size(), std::vector<TGraph *>(3)  );

		// map shower hits to  distance to the closest track hit
		std::vector<std::map< art::Ptr< recob::Hit>, double >> sh_dist(3);
		// vector to save hit with largest minimum distance on each plane
		std::vector< std::pair<art::Ptr< recob::Hit >, double > > max_min_hit(3);

                art::Ptr<recob::PFParticle> pfp_s = showerToPFParticleMap.at(showers[0]);
                auto showerhits = pfParticleToHitsMap.at(pfp_s);

                std::vector<TGraph*> t_pts_s(3);
                std::vector<std::vector<double>>  vec_t(3);
                std::vector<std::vector<double>>  vec_c(3);	
		std::vector<int> num_shr_hits(3);

                for(auto &sh: showerhits){
		    int plane = (int)sh->View();
		    num_shr_hits[plane] += 1;
		    
		    double minDist = 999.9;
		    double dist;
		    // only do if there are track hits on this plane with which to compare
		    if (t_vec_c[(int)sh->View()].size() != 0){	
	                double wire = (double)sh->WireID().Wire;
                        vec_c[(int)sh->View()].push_back(wire); 
			double time = (double)sh->PeakTime();           
			vec_t[(int)sh->View()].push_back(time);

			for (unsigned int th = 0; th < t_vec_c[(int)sh->View()].size(); th++){
			    dist = sqrt( pow (time/25 - (t_vec_t[(int)sh->View()])[th]/25, 2) + pow( wire*0.3 - (t_vec_c[(int)sh->View()])[th]*0.3, 2) );
			    if (dist < minDist) { 
				minDist = dist;
			     }

			} // end of track hits for
			s_hists[(int)sh->View()]->Fill(minDist);

			// keep track of 10 smallest distances and their corresponding hits
			if (sh_dist[plane].size() < 10){
			    (sh_dist[plane])[sh] = minDist;			    // insert shower hit with accompanying minimum distance
			    max_min_hit[plane] = (*std::max_element(sh_dist[plane].begin(), sh_dist[plane].end(), map_max_fn));	    // find new min dist hit
			}
			else{ if (minDist < max_min_hit[plane].second){
				sh_dist[plane].erase(max_min_hit[plane].first);	    // erase old max of min distances from map
				(sh_dist[plane])[sh] = minDist;			    // insert shower hit with accompanying minimum distance	
				max_min_hit[plane] = (*std::max_element(sh_dist[plane].begin(), sh_dist[plane].end(), map_max_fn)); // find new min dist hit
			} }

			// finds the necessary plot boundaries to fit the shower
                        tick_max = std::max(tick_max, (double)sh->PeakTime());
                        tick_min = std::min(tick_min, (double)sh->PeakTime());
                        chan_max[(int)sh->View()] = std::max( chan_max[(int)sh->View()], wire);
                        chan_min[(int)sh->View()] = std::min( chan_min[(int)sh->View()], wire);
		    } // end if stmnt t_vec_c
                } // end looping shower hits

		// create graphs from newly compiled shower coordinates
                t_pts_s[0] = new TGraph(vec_c[0].size(), &(vec_c[0])[0], &(vec_t[0])[0]);
                t_pts_s[1] = new TGraph(vec_c[1].size(), &(vec_c[1])[0], &(vec_t[1])[0]);
                t_pts_s[2] = new TGraph(vec_c[2].size(), &(vec_c[2])[0], &(vec_t[2])[0]);
		// save new graphs for this shower into vector containing all showers
		pts_shr[0] = t_pts_s;
	    
		// place data into approriate vertex_tree variables
		for(int plane = 0; plane < 3; plane++){
		    if (num_shr_hits[plane] == 0){ // if there are no shower hits on this plane, is extremely isolated
			m_isolation_min_dist_trk_shr.push_back(999); 
			m_isolation_nearest_shr_hit_to_trk_wire.push_back(999);
			m_isolation_nearest_shr_hit_to_trk_time.push_back(999);
		    }
		    else if (t_vec_t[plane].size() > 0) { // have to have both shower and track hits on this plane to have valid comparisons for distance
			auto abs_min = (*std::min_element(sh_dist[plane].begin(), sh_dist[plane].end(), map_min_fn));
			m_isolation_min_dist_trk_shr.push_back(abs_min.second); 
			m_isolation_nearest_shr_hit_to_trk_wire.push_back((double)abs_min.first->WireID().Wire);
			m_isolation_nearest_shr_hit_to_trk_time.push_back((double)abs_min.first->PeakTime());
		    }
		    else{ // if there are no shower hits or there are no track hits on this plane, getting min distance fails
			m_isolation_min_dist_trk_shr.push_back(-999); 
			m_isolation_nearest_shr_hit_to_trk_wire.push_back(-999);
			m_isolation_nearest_shr_hit_to_trk_time.push_back(-999);
		    }
		    m_isolation_num_shr_hits_win_1cm_trk.push_back(s_hists[plane]->Integral(1,1));
		    m_isolation_num_shr_hits_win_2cm_trk.push_back(s_hists[plane]->Integral(1,2));
		    m_isolation_num_shr_hits_win_5cm_trk.push_back(s_hists[plane]->Integral(1,5));
		    m_isolation_num_shr_hits_win_10cm_trk.push_back(s_hists[plane]->Integral(1,10));
		}

		/* DRAW SHOWER HISTOGRAM */
		histcan->cd(1);
		s_hists[0]->Draw();
		s_hists[0]->GetXaxis()->SetTitle("distance from shower hit to closest track hit [cm]");
    
		histcan->cd(2);
		s_hists[1]->Draw();
		s_hists[1]->GetXaxis()->SetTitle("distance from shower hit to closest track hit [cm]");

		histcan->cd(3);
		s_hists[2]->Draw();
		s_hists[2]->GetXaxis()->SetTitle("distance from shower hit to closest track hit [cm]");



//NOW THE UNASSOCIATED HITS - that is ones not part of an identified track or shower
std::cout << "Isolation: Acquiring unassociated hits coordinates and comparing with track hits" << std::endl;

	    	// create vector of three layers for unassoc hits 
		std::vector<TGraph *> g_unass(3);

		std::vector<std::vector<std::vector<double>>> pts_to_recluster(3); //plane, point, {wire,tic}
		//make a map tp actual hits here I guess.
		std::vector<std::map<int, art::Ptr<recob::Hit>>> mapPointIndexToHit(3);

		std::vector<double>  minDist_tot(3);
		std::vector<double> minWire(3);
		std::vector<double> minTime(3);

		for(int plane = 0; plane < 3; plane++){
		    minDist_tot[plane] = 999;
		    std::vector<double> vec_t;
		    std::vector<double> vec_c;

		    for(auto &uh: unassociated_hits_all[plane]){
	    
			if (t_vec_c[plane].size() == 0) break;

			double wire = (double)uh->WireID().Wire;
			vec_c.push_back(wire);
			double time = (double)uh->PeakTime();
			vec_t.push_back(time);

			double minDist = 999.9;
			double dist;
			for (unsigned int th = 0; th < t_vec_c[(int)uh->View()].size(); th++){
			    dist = sqrt( pow (time/25 - (t_vec_t[(int)uh->View()])[th]/25, 2) + pow( wire*0.3 - (t_vec_c[(int)uh->View()])[th]*0.3, 2) );
			    if (dist < minDist) { minDist = dist; }
			}
			u_hists[(int)uh->View()]->Fill(minDist);

			if (minDist < minDist_tot[plane]) { 
			    minDist_tot[plane] = minDist;
			    minWire[plane] = wire;
			    minTime[plane] = time;
			}

			// for reclustering
			std::vector<double> pt = {wire, vec_t.back()};
			pts_to_recluster[(int)uh->View()].push_back(pt);
			mapPointIndexToHit[(int)uh->View()][pts_to_recluster[(int)uh->View()].size()-1] = uh;

		    } // end looping unassociated_hits_all

		    g_unass[plane] = new TGraph(vec_c.size(), &vec_c[0], &vec_t[0]);
		} // end looping planes
	
		// place data into appropriate vertex_tree variables
		for(int plane = 0; plane < 3; plane++){
		    if (t_vec_t[plane].size() > 0 && unassociated_hits_all[plane].size() > 0){	
			m_isolation_min_dist_trk_unassoc.push_back(minDist_tot[plane]);  
			m_isolation_nearest_unassoc_hit_to_trk_wire.push_back(minWire[plane]);
			m_isolation_nearest_unassoc_hit_to_trk_time.push_back(minTime[plane]);
		    }
		    else {		
			m_isolation_min_dist_trk_unassoc.push_back(-999);  
			m_isolation_nearest_unassoc_hit_to_trk_wire.push_back(-999);
			m_isolation_nearest_unassoc_hit_to_trk_time.push_back(-999);
		}
		    m_isolation_num_unassoc_hits_win_1cm_trk.push_back(u_hists[plane]->Integral(1,1));
		    m_isolation_num_unassoc_hits_win_2cm_trk.push_back(u_hists[plane]->Integral(1,2));
		    m_isolation_num_unassoc_hits_win_5cm_trk.push_back(u_hists[plane]->Integral(1,5));
		    m_isolation_num_unassoc_hits_win_10cm_trk.push_back(u_hists[plane]->Integral(1,10));		
		}

		/* DRAW UNASSOCIATED HITS DISTANCE HISTOGRAMS */
		histcan->cd(4);
		u_hists[0]->Draw();
		u_hists[0]->GetXaxis()->SetTitle("distance from unassoc hit to closest track hit [cm]");
	    
		histcan->cd(5);
		u_hists[1]->Draw();
		u_hists[1]->GetXaxis()->SetTitle("distance from unassoc hit to closest track hit [cm]");
	    
		histcan->cd(6);
		u_hists[2]->Draw();
		u_hists[2]->GetXaxis()->SetTitle("distance from unassoc hit to closest track hit [cm]");


/* SAVE THE CANVAS -which holds all 6 distance histograms for shower and unassociated hits*/
/*		histcan->Update();
		histcan->SaveAs((print_name + "_hist.pdf").c_str(), "pdf");
*/

		delete histcan;


//PLOTTING NOW
//SET-UP
            double plot_point_size = 0.6;        

            std::vector<int> tcols = {kRed-7, kBlue-7, kGreen-3, kOrange-3, kCyan-3, kMagenta-3, kGreen+1 , kRed+1};
            int used_col=0;
            if(showers.size()+tracks.size() > tcols.size()){
                for(int i =0; i< (int)(showers.size()+tracks.size() - tcols.size() +2); i++){
                    tcols.push_back(tcols[(int)rangen->Uniform(0,7)]+(int)rangen->Uniform(-5,5));
                }
            }

            std::cout<<"*Tick Min: "<<tick_min<<" Max: "<<tick_max<<std::endl;
            auto const TPC = (*geom).begin_TPC();
            auto ID = TPC.ID();
            int fCryostat = ID.Cryostat;
            int fTPC = ID.TPC;
            std::cout<<TPC.ID()<<"*= the beginning TPC ID" <<std::endl;
            std::cout<<"*the cryostat id = "<<fCryostat<<std::endl;  
            std::cout<<"*the tpc id = "<<fTPC<<std::endl;  


//PLOTTING THE VERTEX POSITION ON THE PLOT

            std::vector<double> vertex_time(3); 
            std::vector<double> vertex_wire(3); 


            std::vector<TGraph*> g_vertex(3);
            for(int i=0; i<3; i++){
                TPad * pader = (TPad*)can->cd(i+1);

                if(i==0 ) pader->SetLeftMargin(0.1);

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
                g_vertex[i]->GetYaxis()->SetRangeUser(tick_min*0.9,tick_max*1.1);
                g_vertex[i]->GetXaxis()->SetLimits(chan_min[i]*0.9,chan_max[i]*1.1);
                g_vertex[i]->SetTitle(("Plane " +std::to_string(i)).c_str());
                g_vertex[i]->GetYaxis()->SetTitle("Peak Hit Time Tick");
                g_vertex[i]->GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
                g_vertex[i]->Draw("ap");

                if(i>0){
                    g_vertex[i]->GetYaxis()->SetLabelOffset(999);
                    g_vertex[i]->GetYaxis()->SetLabelSize(0);
                }


            }

// ******************************** DeadWireRegions ********************************************
            for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
                int badchan = bad_channel_list_fixed_mcc9[i].first;                                       
                int ok = bad_channel_list_fixed_mcc9[i].second;       

                if(ok>1)continue;
                auto hs = geom->ChannelToWire(badchan);

                int thisp = (int)hs[0].Plane;
                double bc = hs[0].Wire;
                //                        std::cout<<"WIRE "<<thisp<<" "<<bc<<" "<<hs.size()<<std::endl;
                if(chan_min[thisp]*0.9 < bc && bc < chan_max[thisp]*1.1 ){
                    can->cd(thisp+1);
                    TLine *l = new TLine(bc,tick_min*0.9,bc,tick_max*1.1);
                    l->SetLineColor(kGray+1);
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
		if (g_unass[i]->GetN() > 0){
                    g_unass[i]->SetMarkerColor(kBlack);
                    g_unass[i]->SetMarkerStyle(20);
                    g_unass[i]->SetMarkerSize(plot_point_size);
                }
                g_vertex[i]->Draw("p same");
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

// PLOTTING SHOWER?
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

// PLOTTING 
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

            can->Update();
     //       can->SaveAs((print_name+".pdf").c_str(),"pdf");
            std::cout<<"*PRINTING"<<std::endl;

            delete can;



            }


        }




        return;
    }

}

