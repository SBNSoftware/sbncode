#include "sbncode/SinglePhotonAnalysis/Libraries/second_shower_search.h"

#include "art/Framework/Principal/Event.h"

#include "TPrincipal.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"

#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/Processors.h"

#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_math.h"

namespace single_photon
{

  TGraph* GetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts){

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

  sss_score ScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower){
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
    //std::cout<<"SSS\t||\tEigen: "<<score.pca_0<<" "<<score.pca_1<<std::endl;

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

  int CompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,
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



  std::vector<double>SecondShowerMatching(
      std::vector<art::Ptr<recob::Hit>>& hitz,
      art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
      std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
      //            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
      std::map< int ,art::Ptr<simb::MCParticle>>  & MCParticleToTrackIdMap,
      var_all& vars){


    std::vector<double> ans; //matched,pdg,parentpdg,trkid


    std::vector<double> vec_fraction_matched;
    std::map<std::string,bool> map_is_shower_process = {{"compt",true},{"FastScintillation",true},{"eBrem",true},{"phot",true},{"eIoni",true},{"conv",true},{"annihil",true}};
    bool reco_verbose = false;

    std::unordered_map<int,double> objide; //map between the MCParticle track ID and the backtracker energy

    //energy for an MCParticle that comprises the most energy when sum over associated hits in PFP
    //total energy of the reco PFP taken from the sum of the hits associated to an MCParticle
    // double maxe=-1, tote=0;                // tote unused
    double maxe=-1;                

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
        // tote += match_vec[i_p]->energy; //calculate total energy deposited // unused
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

    //            if(reco_verbose)std::cout << "recoMC()\t||\t On Object "<<i<<". The number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       
    //          if(reco_verbose) std::cout<<"recoMC()\t||\t the fraction of hits from overlay is is "<<fraction_num_hits_overlay<<" ("<<n_not_associated_hits<<"/"<<obj_hits_ptrs.size()<<")"<<std::endl;


    if(n_associated_mcparticle_hits == 0){
      //This will only occur if the whole recob::PFParticle is PURELY associated with an overlay object
      found_a_match =false;
      //Here we will fill every sivars.m_shower_XXX variable with -999 or something like that 
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
    if(reco_verbose)           std::cout<<"recoMC()\t||\t the number of source mother particles is "<<mother_MCP_map.size()<<" of which : "<<marks_mother_vector.size()<<" are unique!"<<std::endl;

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
      //if(reco_verbose)  std::cout<<"recoMC()\t||\t The `BEST` mother is a "<<marks_mother_vector[best_mother_index]->PdgCode()<<" at "<<best_mother_index<<" on plane: "<<best_mother_plane<<std::endl;
      std::cout<<"recoMC()\t||\t The `BEST` mother is a "<<marks_mother_vector[best_mother_index]->PdgCode()<<" at "<<best_mother_index<<" on plane: "<<best_mother_plane<<std::endl;
      for(int l=0; l<3; l++){
        std::cout<<"recoMC()\t||\t It represents "<<marks_mother_energy_fraction_map[marks_mother_vector[best_mother_index]][l]/total_energy_on_plane[l]*100.0<<"% of the energy on plane: "<<l<<" which is "<<total_energy_on_plane[l] <<std::endl;
      }
    }


    if(reco_verbose) std::cout<<"---------------------------- L2-------------------------------"<<std::endl;
    const art::Ptr<simb::MCParticle> match = marks_mother_vector[best_mother_index];

    std::vector<double> corrected_vertex(3), corrected_start(3);
    spacecharge_correction(match, corrected_vertex);


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







  void SecondShowerSearch3D(
      std::vector<art::Ptr<recob::Shower>> & showers,
      std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & NormalShowerToPFParticleMap,  
      std::vector<art::Ptr<recob::Track>> & tracks, 
      std::map<art::Ptr<recob::Track>,  
      art::Ptr<recob::PFParticle>> & NormalTrackToPFParticleMap, 
      art::Event const & evt, 
      var_all& vars,
      para_all& paras
      ){

    std::string sss3dlabel = "pandoraShower";//"pandoraAllOutcomesShower" WARNING, should check this!

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
    vars.m_sss3d_num_showers = (int)n_all_shr-showers.size();

    if(showers.size()==0) return;

    auto primary_shower = showers.front();

    std::cout<<"PandoraAllOutcomesShower has "<<n_all_shr<<" Showers in total. "<<std::endl;
    for(auto &shr: allShowerVector){
      //lets look at 3D distance to "vertex", and only go further with things that are within 80cm [default]
      double dist = sqrt(pow(vars.m_vertex_pos_x - shr->ShowerStart().X(),2)+pow(vars.m_vertex_pos_y - shr->ShowerStart().Y(),2)+pow(vars.m_vertex_pos_z - shr->ShowerStart().Z(),2) );
      if(dist>paras.s_max_conv_dist) continue;

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

      //std::cout<<shr.key()<<" Length "<<shr->Length()<<" Dist "<<dist<<" Impact: "<<impact_paramater_shr(vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z, shr)<<" pfp key: "<<showerToPFParticleMap[shr].key()<<" Self "<<showerToPFParticleMap[shr]->Self()<<std::endl;

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

      double senergy = std::max( CalcEShowerPlane(showerToHitsMap[shr], 0, paras), std::max( CalcEShowerPlane(showerToHitsMap[shr], 1, paras),  CalcEShowerPlane(showerToHitsMap[shr], 2, paras)) );
      double invar = implied_invar_mass(vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z,  primary_shower, vars.m_reco_shower_energy_max[0], shr, senergy);
      double implied_invar = invar_mass(primary_shower, vars.m_reco_shower_energy_max[0], shr, senergy) ;
      double shr_score = 0.0; //need pfp and metadata to get score, and might give slice! (This will be harder..) but on reflection, kinda important. PCA spread might be a good rplacement.
      int is_clear_cosmic_slice = 0 ;
      int is_nu_slice = 0;


      vars.m_sss3d_shower_start_x.push_back(shr->ShowerStart().X());
      vars.m_sss3d_shower_start_y.push_back(shr->ShowerStart().Y());
      vars.m_sss3d_shower_start_z.push_back(shr->ShowerStart().Z());
      vars.m_sss3d_shower_dir_x.push_back(shr->Direction().X());
      vars.m_sss3d_shower_dir_y.push_back(shr->Direction().Y());
      vars.m_sss3d_shower_dir_z.push_back(shr->Direction().Z());
      vars.m_sss3d_shower_length.push_back(shr->Length());
      vars.m_sss3d_shower_conversion_dist.push_back(dist);
      vars.m_sss3d_shower_invariant_mass.push_back(invar);
      vars.m_sss3d_shower_implied_invariant_mass.push_back(implied_invar);
      double imp = impact_paramater_shr(vars.m_vertex_pos_x, vars.m_vertex_pos_y, vars.m_vertex_pos_z, shr);
      vars.m_sss3d_shower_impact_parameter.push_back(imp);

      if(dist!=0) {
        vars.m_sss3d_shower_ioc_ratio.push_back(imp/dist);
      }else{
        vars.m_sss3d_shower_ioc_ratio.push_back(0);

      }
      vars.m_sss3d_shower_energy_max.push_back(senergy);
      vars.m_sss3d_shower_score.push_back(shr_score);
      vars.m_sss3d_slice_clear_cosmic.push_back(is_clear_cosmic_slice);
      vars.m_sss3d_slice_nu.push_back(is_nu_slice);
    }

    return;
  }



  void SimpleSecondShowerCluster(var_all& vars, para_all& paras){

    std::string base = "sss3d_";
    std::vector<std::string> mod = {"ioc_ranked","invar_ranked"};

    vars.m_sss3d_ioc_ranked_en = -9;
    vars.m_sss3d_ioc_ranked_conv = -9;
    vars.m_sss3d_ioc_ranked_invar = -9;
    vars.m_sss3d_ioc_ranked_implied_invar = -9;
    vars.m_sss3d_ioc_ranked_ioc = -9;
    vars.m_sss3d_ioc_ranked_opang = -9;
    vars.m_sss3d_ioc_ranked_implied_opang = -9; 
    vars.m_sss3d_ioc_ranked_id = -9;

    vars.m_sss3d_invar_ranked_en = -9;
    vars.m_sss3d_invar_ranked_conv = -9;
    vars.m_sss3d_invar_ranked_invar = -9;
    vars.m_sss3d_invar_ranked_implied_invar = -9;
    vars.m_sss3d_invar_ranked_ioc = -9;
    vars.m_sss3d_invar_ranked_opang = -9;
    vars.m_sss3d_invar_ranked_implied_opang = -9; 
    vars.m_sss3d_invar_ranked_id = -9;


    std::string base2d = "sss_";
    std::vector<std::string> mod2d = {"ioc_ranked","conv_ranked","invar_ranked"};

    vars.m_sss2d_ioc_ranked_en = -9;
    vars.m_sss2d_ioc_ranked_conv = -9;
    vars.m_sss2d_ioc_ranked_ioc = -9;
    vars.m_sss2d_ioc_ranked_pca = -9;
    vars.m_sss2d_ioc_ranked_invar = -9;
    vars.m_sss2d_ioc_ranked_angle_to_shower = -9;
    vars.m_sss2d_ioc_ranked_num_planes = -9;

    vars.m_sss2d_conv_ranked_en = -9;
    vars.m_sss2d_conv_ranked_conv = -9;
    vars.m_sss2d_conv_ranked_ioc = -9;
    vars.m_sss2d_conv_ranked_pca = -9;
    vars.m_sss2d_conv_ranked_invar = -9;
    vars.m_sss2d_conv_ranked_angle_to_shower = -9;
    vars.m_sss2d_conv_ranked_num_planes = -9;

    vars.m_sss2d_invar_ranked_en = -9;
    vars.m_sss2d_invar_ranked_conv = -9;
    vars.m_sss2d_invar_ranked_ioc = -9;
    vars.m_sss2d_invar_ranked_pca = -9;
    vars.m_sss2d_invar_ranked_invar = -9;
    vars.m_sss2d_invar_ranked_angle_to_shower = -9;
    vars.m_sss2d_invar_ranked_num_planes = -9;

    //---------------------------------------
    //First off, the 3D showers 
    //First some 3D shower information
    if(vars.m_sss3d_shower_conversion_dist.size()>0 && vars.m_reco_shower_energy_max.size()>0){
      //std::cout<<"Primary shower en "<<reco_shower_energy_max->at(0)<<std::endl;

      std::vector<double> inv = vars.m_sss3d_shower_implied_invariant_mass;
      for(auto &v : inv) v = fabs(v-paras.s_mass_pi0_mev);

      std::vector<size_t> ranked_ioc = sort_indexes_rev<double>((vars.m_sss3d_shower_ioc_ratio));
      std::vector<size_t> ranked_invar = sort_indexes_rev<double>((inv));
      std::vector<size_t> ranked_conv = sort_indexes_rev<double>((vars.m_sss3d_shower_conversion_dist));
      std::vector<size_t> ranked_en = sort_indexes_rev<double>((vars.m_sss3d_shower_energy_max));

      int to_consider = vars.m_sss3d_shower_conversion_dist.size();

      if(false){
        std::cout<<"IOC"<<std::endl;
        for(int j=0; j<to_consider; j++){
          std::cout<<"--"<<ranked_ioc[j]<<" ioc: "<<vars.m_sss3d_shower_ioc_ratio.at( ranked_ioc[j] )<<" invar: "<<vars.m_sss3d_shower_implied_invariant_mass.at(ranked_ioc[j])<<" en: "<<vars.m_sss3d_shower_energy_max.at(ranked_ioc[j])<<" conv: "<<vars.m_sss3d_shower_conversion_dist.at(ranked_ioc[j])<<std::endl;
        }

        std::cout<<"INVAR"<<std::endl;
        for(int j=0; j<to_consider; j++){
          std::cout<<"--"<<ranked_invar[j]<<" ioc: "<<vars.m_sss3d_shower_ioc_ratio.at( ranked_invar[j] )<<" invar: "<<vars.m_sss3d_shower_implied_invariant_mass.at(ranked_invar[j])<<" en: "<<vars.m_sss3d_shower_energy_max.at(ranked_invar[j])<<" conv: "<<vars.m_sss3d_shower_conversion_dist.at(ranked_invar[j]) <<std::endl;
        }

        std::cout<<"EN"<<std::endl;
        for(int j=0; j<to_consider; j++){
          int rk = ranked_en[ranked_en.size()-1-j];
          std::cout<<"--"<<rk<<" ioc: "<<vars.m_sss3d_shower_ioc_ratio.at( rk )<<" invar: "<<vars.m_sss3d_shower_implied_invariant_mass.at(rk)<<" en: "<<vars.m_sss3d_shower_energy_max.at(rk)<<" conv: "<<vars.m_sss3d_shower_conversion_dist.at(rk)<<std::endl;
        }
        std::cout<<"CONV"<<std::endl;
        for(int j=0; j<to_consider; j++){
          std::cout<<"--"<<ranked_conv[j]<<" ioc: "<<vars.m_sss3d_shower_ioc_ratio.at( ranked_conv[j] )<<" invar: "<<vars.m_sss3d_shower_implied_invariant_mass.at(ranked_conv[j])<<" en: "<<vars.m_sss3d_shower_energy_max.at(ranked_conv[j])<<" conv: "<<vars.m_sss3d_shower_conversion_dist.at(ranked_conv[j])<<std::endl;
        }
      }


      //std::cout<<"Best IOC "<<"--"<<ranked_ioc[0]<<" ioc: "<<sss3d_shower_ioc_ratio->at( ranked_ioc[0] )<<" invar: "<<sss3d_shower_implied_invariant_mass->at(ranked_ioc[0])<<" en: "<<sss3d_shower_energy_max->at(ranked_ioc[0])<<" conv: "<<sss3d_shower_conversion_dist->at(ranked_ioc[0])<<std::endl;
      vars.m_sss3d_ioc_ranked_en = vars.m_sss3d_shower_energy_max.at(ranked_ioc[0]);            
      vars.m_sss3d_ioc_ranked_conv = vars.m_sss3d_shower_conversion_dist.at(ranked_ioc[0]);            
      vars.m_sss3d_ioc_ranked_invar = vars.m_sss3d_shower_invariant_mass.at(ranked_ioc[0]);            
      vars.m_sss3d_ioc_ranked_implied_invar = vars.m_sss3d_shower_implied_invariant_mass.at(ranked_ioc[0]);            
      vars.m_sss3d_ioc_ranked_ioc = vars.m_sss3d_shower_ioc_ratio.at(ranked_ioc[0]);
      vars.m_sss3d_ioc_ranked_opang = 1.0 - pow(vars.m_sss3d_shower_invariant_mass.at(ranked_ioc[0]),2)/(2.0*vars.m_sss3d_shower_energy_max.at(ranked_ioc[0])*vars.m_reco_shower_energy_max.at(0));
      vars.m_sss3d_ioc_ranked_implied_opang = 1.0 - pow(vars.m_sss3d_shower_implied_invariant_mass.at(ranked_ioc[0]),2)/(2.0*vars.m_sss3d_shower_energy_max.at(ranked_ioc[0])*vars.m_reco_shower_energy_max.at(0));
      vars.m_sss3d_ioc_ranked_id = ranked_ioc[0];


      // std::cout<<"Best invar "<<"--"<<ranked_invar[0]<<" ioc: "<<sss3d_shower_ioc_ratio.at( ranked_invar[0] )<<" invar: "<<sss3d_shower_implied_invariant_mass.at(ranked_invar[0])<<" en: "<<sss3d_shower_energy_max.at(ranked_invar[0])<<" conv: "<<sss3d_shower_conversion_dist.at(ranked_invar[0])<<std::endl;

      // minimum discrepancy between implied invariant mass and pi0 mass
      vars.m_sss3d_invar_ranked_en = vars.m_sss3d_shower_energy_max.at(ranked_invar[0]);            
      vars.m_sss3d_invar_ranked_conv = vars.m_sss3d_shower_conversion_dist.at(ranked_invar[0]);            
      vars.m_sss3d_invar_ranked_invar = vars.m_sss3d_shower_invariant_mass.at(ranked_invar[0]);            
      vars.m_sss3d_invar_ranked_implied_invar = vars.m_sss3d_shower_implied_invariant_mass.at(ranked_invar[0]);            
      vars.m_sss3d_invar_ranked_ioc = vars.m_sss3d_shower_ioc_ratio.at(ranked_invar[0]);
      vars.m_sss3d_invar_ranked_opang = 1.0 - pow(vars.m_sss3d_shower_invariant_mass.at(ranked_invar[0]),2)/(2.0*vars.m_sss3d_shower_energy_max.at(ranked_invar[0])*vars.m_reco_shower_energy_max.at(0));
      vars.m_sss3d_invar_ranked_implied_opang = 1.0 - pow(vars.m_sss3d_shower_implied_invariant_mass.at(ranked_invar[0]),2)/(2.0*vars.m_sss3d_shower_energy_max.at(ranked_invar[0])*vars.m_reco_shower_energy_max.at(0));
      vars.m_sss3d_invar_ranked_id = ranked_invar[0];

    }//end of 3D shower searching


    //Now some 2D shower information
    //
    if(vars.m_sss_num_candidates>0){
      //std::cout<<"2D clusters: "<<sss_num_candidates<<std::endl;
      std::vector<int> nplans(3,0);
      std::vector<std::vector<int>> indexmap(3);


      for(int i=0; i< vars.m_sss_num_candidates; i++){
        //std::cout<<i<<" p: "<<vars.m_sss_candidate_plane->at(i)<<" pdg: "<<vars.m_sss_candidate_parent_pdg->at(i)<<" ovf "<<vars.m_sss_candidate_overlay_fraction->at(i)<<" conv: "<<vars.m_sss_candidate_min_dist->at(i)<<std::endl;
      }

      std::vector<std::vector<int>> uniq_candidates;

      for(int i=0; i< vars.m_sss_num_candidates; i++){
        int ip = vars.m_sss_candidate_plane.at(i);
        //int nhits = sss_candidate_num_hits.at(i);
        nplans[ip]++;
        indexmap[ip].push_back(i);

        //Two passes to build up all "Candidates" for 2 and 3 plane matches
        for(int j=i;j<vars.m_sss_num_candidates;j++){
          int jp = vars.m_sss_candidate_plane.at(j);
          if(jp==ip) continue;

          bool contain_ij = false;
          bool contain_ji = false;
          if(vars.m_sss_candidate_mean_tick.at(j)<=vars.m_sss_candidate_max_tick.at(i) && vars.m_sss_candidate_mean_tick.at(j) >= vars.m_sss_candidate_min_tick.at(i))contain_ij = true;
          if(vars.m_sss_candidate_mean_tick.at(i)<=vars.m_sss_candidate_max_tick.at(j) && vars.m_sss_candidate_mean_tick.at(i) >= vars.m_sss_candidate_min_tick.at(j))contain_ji = true;
          //                                std::cout<<i<<" "<<j<<" "<<contain_ij<<" "<<contain_ji<<std::endl;
          if(contain_ij && contain_ji){
            uniq_candidates.push_back({i,j});
          }
        }
      }

      //Now loop over to check if any indlude a third plane
      for(int i = 0; i< (int)uniq_candidates.size(); i++){
        for(int k=0; k<vars.m_sss_num_candidates; k++){
          //first check if this possible 3rd match is on a seperate plane
          bool new_plane = true;
          for(auto &pp:uniq_candidates[i]){
            if(vars.m_sss_candidate_plane.at(k)==vars.m_sss_candidate_plane.at(pp) ) new_plane = false;   
          }
          if(new_plane){

            bool contain_ik = false;
            bool contain_ki = false;
            bool contain_jk = false;
            bool contain_kj = false;
            if(vars.m_sss_candidate_mean_tick.at(k)<=vars.m_sss_candidate_max_tick.at(uniq_candidates[i][0]) && vars.m_sss_candidate_mean_tick.at(k) >= vars.m_sss_candidate_min_tick.at(uniq_candidates[i][0]))contain_ik = true;
            if(vars.m_sss_candidate_mean_tick.at(uniq_candidates[i][0])<=vars.m_sss_candidate_max_tick.at(k) && vars.m_sss_candidate_mean_tick.at(uniq_candidates[i][0]) >= vars.m_sss_candidate_min_tick.at(k))contain_ki = true;
            if(vars.m_sss_candidate_mean_tick.at(k)<=vars.m_sss_candidate_max_tick.at(uniq_candidates[i][1]) && vars.m_sss_candidate_mean_tick.at(k) >= vars.m_sss_candidate_min_tick.at(uniq_candidates[i][1]))contain_ik = true;
            if(vars.m_sss_candidate_mean_tick.at(uniq_candidates[i][1])<=vars.m_sss_candidate_max_tick.at(k) && vars.m_sss_candidate_mean_tick.at(uniq_candidates[i][1]) >= vars.m_sss_candidate_min_tick.at(k))contain_ki = true;

            //If this matches well with Either last candidate, include as a possibility
            if((contain_ik&&contain_ki) || (contain_jk&&contain_kj)){
              uniq_candidates[i].push_back(k);
            }

          }
        }
      }
      //Check which candidates have been used where
      std::vector<int> used_candidates(vars.m_sss_num_candidates);
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
      std::cout<<"After: used_candidates "<<vars.m_sss_num_candidates<<std::endl;
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
        // double mean_conv = 0.0; // unused
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

          //std::cout<<"----- plane: "<<vars.m_sss_candidate_plane.at(ic)<<" nhits: "<<vars.m_sss_candidate_num_hits.at(ic)<<" imp: "<<vars.m_sss_candidate_impact_parameter.at(ic)<<" en: "<<vars.m_sss_candidate_energy.at(ic)<<" a2s: "<<vars.m_sss_candidate_angle_to_shower.at(ic)<<" conv: "<<vars.m_sss_candidate_min_dist.at(ic)<<" pdg: "<<vars.m_sss_candidate_parent_pdg.at(ic)<<std::endl;

          double eff_invar = sqrt(2.0*vars.m_sss_candidate_energy.at(ic)*vars.m_reco_shower_energy_max.at(0)*(1.0-cos(vars.m_sss_candidate_angle_to_shower.at(ic))));
          double eff_invar_diff = fabs(eff_invar-139.5);

          //is_min_slice = std::min(is_in_slice, vars.m_sss_candidate_in_nu_slice[ic] );
          //is_in_slice.push_back(vars.m_sss_candidate_in_nu_slice[ic]); 

          mean_min_dist +=vars.m_sss_candidate_min_dist.at(ic)/(double)nt;
          mean_energy +=vars.m_sss_candidate_energy.at(ic)/(double)nt;
          mean_impact +=vars.m_sss_candidate_impact_parameter.at(ic)/(double)nt;
          // mean_conv +=vars.m_sss_candidate_min_dist.at(ic)/(double)nt; // unused
          mean_invar +=eff_invar/(double)nt;
          mean_invar_diff +=eff_invar_diff/(double)nt;

          min_conv = std::min(min_conv, vars.m_sss_candidate_min_dist.at(ic));
          max_min_dist = std::max(max_min_dist, vars.m_sss_candidate_min_dist.at(ic));
          max_pca = std::max(max_pca, vars.m_sss_candidate_PCA.at(ic));
          min_impact = std::min(min_impact, vars.m_sss_candidate_impact_parameter.at(ic));
          min_angle = std::min(min_angle, vars.m_sss_candidate_angle_to_shower.at(ic));

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

      vars.m_sss2d_ioc_ranked_en = candidates_en[ranked_ioc[0]];
      vars.m_sss2d_ioc_ranked_conv = candidates_conv[ranked_ioc[0]];
      vars.m_sss2d_ioc_ranked_ioc = candidates_ioc[ranked_ioc[0]];
      vars.m_sss2d_ioc_ranked_invar = candidates_eff_invar[ranked_ioc[0]];
      vars.m_sss2d_ioc_ranked_pca = candidates_pca[ranked_ioc[0]];
      vars.m_sss2d_ioc_ranked_angle_to_shower  = candidates_angle_to_shower[ranked_ioc[0]];
      vars.m_sss2d_ioc_ranked_num_planes  = candidates_num_planes[ranked_ioc[0]];

      vars.m_sss2d_conv_ranked_en = candidates_en[ranked_conv[0]];
      vars.m_sss2d_conv_ranked_conv = candidates_conv[ranked_conv[0]];
      vars.m_sss2d_conv_ranked_ioc = candidates_ioc[ranked_conv[0]];
      vars.m_sss2d_conv_ranked_invar = candidates_eff_invar[ranked_conv[0]];
      vars.m_sss2d_conv_ranked_pca = candidates_pca[ranked_conv[0]];
      vars.m_sss2d_conv_ranked_angle_to_shower  = candidates_angle_to_shower[ranked_conv[0]];
      vars.m_sss2d_conv_ranked_num_planes  = candidates_num_planes[ranked_conv[0]];

      vars.m_sss2d_invar_ranked_en = candidates_en[ranked_invar[0]];
      vars.m_sss2d_invar_ranked_conv = candidates_conv[ranked_invar[0]];
      vars.m_sss2d_invar_ranked_ioc = candidates_ioc[ranked_invar[0]];
      vars.m_sss2d_invar_ranked_invar = candidates_eff_invar[ranked_invar[0]];
      vars.m_sss2d_invar_ranked_pca = candidates_pca[ranked_invar[0]];
      vars.m_sss2d_invar_ranked_angle_to_shower  = candidates_angle_to_shower[ranked_invar[0]];
      vars.m_sss2d_invar_ranked_num_planes  = candidates_num_planes[ranked_invar[0]];
    }

    return ;
  }



  std::pair<bool, std::vector<double>> clusterCandidateOverlap(const std::vector<int> & candidate_indices, const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks){

    size_t size = candidate_indices.size();
    if(size == 0){
      throw std::runtime_error("clusterCandidateOverlap: No cluster candidates to analyze time overlap for..");
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


  std::pair<int, std::pair<std::vector<std::vector<double>>, std::vector<double>>> GroupClusterCandidate(int num_clusters,  const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks){
    std::cout << "group_cluster_candidate\t|| Total of " << num_clusters << " to be grouped" << std::endl;

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

    std::cout << "GroupClusterCandidate\t|| Formed " << num_cluster_groups << " cluster groups" << std::endl;

    return {num_cluster_groups, {grouped_cluster_indices, cluster_group_timeoverlap_fraction}};
  } 

  //isolation.h
  /* Arguments to Function  IsolationStudy  (all are const):
   * 1. vector named tracks     of art ptr to recob track
   * 2. map named trackToPFPParticleMap   of .i. art ptr to recob track    .ii. art ptr to recob pfparticle  
   * 3. vector named showers           of art ptr to recob showers
   * 4. map named showerToPFParticleMap   of .i. art ptr to recob showe    .ii. art ptr to recob pfparticle
   * 5. map named pfParticleToHistMap     of .i. art ptr to recob prparticle  .ii. vector of art ptr to recob hit
   * 6. map named pfParticleToSliceIDMap  of .i. art ptr to recob pfparticle  .ii. int
   * 7. map named sliceIDToHitsMap  of .i. int and        .ii. art ptr to recob hit
   */
  void IsolationStudy(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Track>>& tracks, 
      const std::vector<art::Ptr<recob::Shower>>& showers, 
      detinfo::DetectorPropertiesData const & theDetector,
      var_all& vars,
      para_all& paras){

    int total_track_hits =0;
    int total_shower_hits =0;
    int nu_slice_id = -999;

    std::vector< art::Ptr<recob::Hit> > associated_hits;
    std::vector< art::Ptr<recob::Hit> > unassociated_hits;
    std::vector< art::Ptr<recob::Hit> > unassociated_hits_plane0;
    std::vector< art::Ptr<recob::Hit> > unassociated_hits_plane1;
    std::vector< art::Ptr<recob::Hit> > unassociated_hits_plane2;

    std::vector< std::map<size_t, std::vector< art::Ptr<recob::Hit> >> > v_newClusterToHitsMap(3);

    std::vector<art::Ptr<recob::Hit>> slicehits;

    // BEGIN FOR LOOP TO COUNT TRACK HITS AND PUSH INTO ASSOCIATED HITS
    for(size_t t =0; t< tracks.size(); t++){
      art::Ptr<recob::Track> track = tracks[t];
      PandoraPFParticle* ppfp = PPFP_GetPPFPFromTrack(all_PPFPs, track);
      //      art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;//trackToPFParticleMap[track];

      int sliceid = ppfp->get_SliceID();//pfParticleToSliceIDMap.at(pfp);
      //WARNING the temp. solution only work with best nuscore slice Keng
      if(!ppfp->get_IsNuSlice()) continue;

      std::vector<art::Ptr<recob::Hit>> tmp_slicehits = ppfp->pSliceHits;//sliceIDToHitsMap.at(sliceid);
      std::vector<art::Ptr<recob::Hit>> trackhits = ppfp->pPFPHits;//pfParticleToHitsMap.at(pfp);

      if(ppfp->get_IsNeutrino()) slicehits.insert(slicehits.end(), tmp_slicehits.begin(), tmp_slicehits.end());//add up all nu slice hits

      std::cout << "*SSS: track "<< t <<" is in slice "<< sliceid <<" which has "<< tmp_slicehits.size() <<" hits. This track has  "<< trackhits.size() <<" of them. " << std::endl;
      total_track_hits += trackhits.size();

      //WARNING  nu_slice_id is updated? Oh, tracks[0] and tracks[1] have different slide IDs (both nu slice but different nu score), need to fix this;
      //       temporary solution is to skip tracks[1];
      //      if(nu_slice_id !=  sliceid && nu_slice_id != -999){
      //        std::cout<<"*ERROR!! In Second Shower Search (Isolation Study), the neutrino slice ID changed? this: "<<sliceid<<", last: "<<nu_slice_id<<std::endl;
      //        exit(EXIT_FAILURE);
      //      }   
      //      nu_slice_id = sliceid;

      for(auto &h: trackhits){
        associated_hits.push_back(h);
      }
    }
    // END FOR LOOPING TRACKS



    // BEGIN FOR LOOPING SHOWERS TO COUNT SHOWER HITS AND PUSH INTO ASSOC HITS
    for(size_t s =0; s< showers.size(); s++){
      art::Ptr<recob::Shower> shower = showers[s];
      PandoraPFParticle* ppfp = PPFP_GetPPFPFromShower(all_PPFPs, shower);
      //            art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;//showerToPFParticleMap.at(shower);

      int sliceid = ppfp->get_SliceID();//pfParticleToSliceIDMap.at(pfp);

      //    if(sliceid != slice_w_bestnuID) continue;//WARNING only deal with nu slice with best nu score for now Keng
      if(!ppfp->get_IsNuSlice()) continue;
      if(sliceid<0) continue; //negative sliceid is bad

      auto tmp_slicehits = ppfp->pSliceHits;//sliceIDToHitsMap.at(sliceid);
      auto showerhits = ppfp->pPFPHits;//pfParticleToHitsMap.at(pfp);

      std::cout<<"*SSS: shower "<< s <<" is in slice "<< sliceid <<" which has "<< tmp_slicehits.size() <<" hits. This shower has  "<< showerhits.size() <<" of them. "<< std::endl;
      total_shower_hits+=showerhits.size();

      for(auto &h: showerhits){
        associated_hits.push_back(h);
      } 
    } 
    // END FOR LOOP COUNTING SHOWER HITS

    vars.m_sss_num_associated_hits = total_shower_hits + total_track_hits;

    // PRINT SUMMARY OF HIT TYPES
    std::cout<<"*SSS: So in total we have "<<total_shower_hits<<" shower hits and "<<total_track_hits<<" track hits"<<" associatedHits total: "<<associated_hits.size()<<std::endl;


    // IF VALID SLICE
    if(nu_slice_id >= 0){
      vars.m_sss_num_unassociated_hits = slicehits.size()-total_shower_hits-total_track_hits;

      std::cout<<"*SSS: So that leaves "<<slicehits.size()-total_shower_hits-total_track_hits<<" hits not included in tracks and showers"<<std::endl;

      if(vars.m_sss_num_unassociated_hits < 0){
        std::cout<<"ERROR!! Number of unassociated hits is negative, i.e: num_associated: "<<vars.m_sss_num_associated_hits<<" and total slice hits: "<<slicehits.size()<<std::endl;
        exit(EXIT_FAILURE);
      }



      // DETERMINE UNASSOCIATED HITS BY COMPARING ALL SLICE HITS WITH LIST OF ASSOCIATED HITS
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
      if(paras.s_make_sss_plots && showers.size()==1 && tracks.size()==1){

        std::string print_name = "isolation_"+std::to_string(vars.m_run_number)+"_"+std::to_string(vars.m_subrun_number)+"_"+std::to_string(vars.m_event_number);
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

        PandoraPFParticle* ppfpt = PPFP_GetPPFPFromTrack(all_PPFPs, tracks[0]);
        //        art::Ptr<recob::PFParticle> pfpt = trackToPFParticleMap.at(tracks[0]);
        auto trackhits = ppfpt->pPFPHits;//pfParticleToHitsMap.at(pfpt);

        std::vector<TGraph*> t_pts(3);
        std::vector<std::vector<double>>  t_vec_t(3);  // peak time of track hits on 3 planes.
        std::vector<std::vector<double>>  t_vec_c(3); // wire number of track hits on 3 planes.

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

        // map shower hits to  distance to the closest track hit (hold up to 10 hits with minimum distance)
        std::vector<std::map< art::Ptr< recob::Hit>, double >> sh_dist(3);
        // vector to save hit with largest minimum distance (in sh_dist) on each plane
        std::vector< std::pair<art::Ptr< recob::Hit >, double > > max_min_hit(3);

        PandoraPFParticle* ppfps = PPFP_GetPPFPFromShower(all_PPFPs, showers[0]);
        auto showerhits = ppfps->pPFPHits;//pfParticleToHitsMap.at(pfp_s);
        //        art::Ptr<recob::PFParticle> pfp_s = showerToPFParticleMap.at(showers[0]);

        std::vector<TGraph*> t_pts_s(3);
        std::vector<std::vector<double>>  vec_t(3);
        std::vector<std::vector<double>>  vec_c(3);  
        std::vector<int> num_shr_hits(3);

        for(auto &sh: showerhits){
          int plane = (int)sh->View();
          num_shr_hits[plane] += 1;

          double minDist = 999.9;  //minimum distance between this shower hit and all track hits 
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
              (sh_dist[plane])[sh] = minDist;          // insert shower hit with accompanying minimum distance
              max_min_hit[plane] = (*std::max_element(sh_dist[plane].begin(), sh_dist[plane].end(), map_max_fn));      // find new min dist hit
            }
            else{ if (minDist < max_min_hit[plane].second){
              sh_dist[plane].erase(max_min_hit[plane].first);      // erase old max of min distances from map
              (sh_dist[plane])[sh] = minDist;          // insert shower hit with accompanying minimum distance  
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
            vars.m_isolation_min_dist_trk_shr.push_back(999); 
            vars.m_isolation_nearest_shr_hit_to_trk_wire.push_back(999);
            vars.m_isolation_nearest_shr_hit_to_trk_time.push_back(999);
          }
          else if (t_vec_t[plane].size() > 0) { // have to have both shower and track hits on this plane to have valid comparisons for distance
            auto abs_min = (*std::min_element(sh_dist[plane].begin(), sh_dist[plane].end(), map_min_fn));
            vars.m_isolation_min_dist_trk_shr.push_back(abs_min.second); 
            vars.m_isolation_nearest_shr_hit_to_trk_wire.push_back((double)abs_min.first->WireID().Wire);
            vars.m_isolation_nearest_shr_hit_to_trk_time.push_back((double)abs_min.first->PeakTime());
          }
          else{ // if there are no shower hits or there are no track hits on this plane, getting min distance fails
            vars.m_isolation_min_dist_trk_shr.push_back(-999); 
            vars.m_isolation_nearest_shr_hit_to_trk_wire.push_back(-999);
            vars.m_isolation_nearest_shr_hit_to_trk_time.push_back(-999);
          }
          vars.m_isolation_num_shr_hits_win_1cm_trk.push_back(s_hists[plane]->Integral(1,1));
          vars.m_isolation_num_shr_hits_win_2cm_trk.push_back(s_hists[plane]->Integral(1,2));
          vars.m_isolation_num_shr_hits_win_5cm_trk.push_back(s_hists[plane]->Integral(1,5));
          vars.m_isolation_num_shr_hits_win_10cm_trk.push_back(s_hists[plane]->Integral(1,10));
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

          if (t_vec_c[plane].size() == 0) break; // if track doesn't have hits on that plane

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

          if (minDist < minDist_tot[plane]) { // of all unassociated hits, the one that has min distance to closest track hits 
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
          vars.m_isolation_min_dist_trk_unassoc.push_back(minDist_tot[plane]);  
          vars.m_isolation_nearest_unassoc_hit_to_trk_wire.push_back(minWire[plane]);
          vars.m_isolation_nearest_unassoc_hit_to_trk_time.push_back(minTime[plane]);
        }
        else {    
          vars.m_isolation_min_dist_trk_unassoc.push_back(-999);  
          vars.m_isolation_nearest_unassoc_hit_to_trk_wire.push_back(-999);
          vars.m_isolation_nearest_unassoc_hit_to_trk_time.push_back(-999);
        }
        vars.m_isolation_num_unassoc_hits_win_1cm_trk.push_back(u_hists[plane]->Integral(1,1));
        vars.m_isolation_num_unassoc_hits_win_2cm_trk.push_back(u_hists[plane]->Integral(1,2));
        vars.m_isolation_num_unassoc_hits_win_5cm_trk.push_back(u_hists[plane]->Integral(1,5));
        vars.m_isolation_num_unassoc_hits_win_10cm_trk.push_back(u_hists[plane]->Integral(1,10));    
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
      /*    histcan->Update();
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
          tcols.push_back(tcols[(int)paras.rangen->Uniform(0,7)]+(int)paras.rangen->Uniform(-5,5));
        }
      }

      std::cout<<"*Tick Min: "<<tick_min<<" Max: "<<tick_max<<std::endl;
      auto const& TPC = *paras.s_geom->begin<geo::TPCGeo>();
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

        std::vector<double> wire = {(double)calcWire(vars.m_vertex_pos_y, vars.m_vertex_pos_z, i, fTPC, fCryostat, *paras.s_geom)};
        std::vector<double> time = {theDetector.ConvertXToTicks(vars.m_vertex_pos_x, i, fTPC,fCryostat)};

        vertex_time[i] = time[0];
        vertex_wire[i] = wire[0];

        if(i==0) vars.m_vertex_pos_wire_p0 = wire[0];
        if(i==1) vars.m_vertex_pos_wire_p1 = wire[0];
        if(i==2) vars.m_vertex_pos_wire_p2 = wire[0];
        vars.m_vertex_pos_tick = time[0];

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
      //plot dead wire  
      //      for(size_t i=0; i< bad_channel_list_fixed_mcc9.size(); i++){
      //        int badchan = bad_channel_list_fixed_mcc9[i].first;                                       
      //        int ok = bad_channel_list_fixed_mcc9[i].second;       
      //
      //        if(ok>1)continue;
      //        auto hs = geom->ChannelToWire(badchan);
      //
      //        int thisp = (int)hs[0].Plane;
      //        double bc = hs[0].Wire;
      //        //                        std::cout<<"WIRE "<<thisp<<" "<<bc<<" "<<hs.size()<<std::endl;
      //        if(chan_min[thisp]*0.9 < bc && bc < chan_max[thisp]*1.1 ){
      //          can->cd(thisp+1);
      //          TLine *l = new TLine(bc,tick_min*0.9,bc,tick_max*1.1);
      //          l->SetLineColor(kGray+1);
      //          l->Draw("same");
      //        }
      //      }

      // plot track
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

      // plot shower hits
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


      // plot unassociated hits
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
      std::string pot_draw = "Run: "+std::to_string(vars.m_run_number)+" SubRun: "+std::to_string(vars.m_subrun_number)+" Event: "+std::to_string(vars.m_event_number);
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
