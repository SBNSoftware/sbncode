#include "sbncode/SinglePhotonAnalysis/Libraries/reco_truth_matching.h"

#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/Processors.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_gadget.h"

namespace single_photon
{
  void RecoMCTracks(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Track>>& tracks,  
      std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,
      std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,
      std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector,  
      std::map< int, art::Ptr<simb::MCParticle> > &      MCParticleToTrackIdMap, 
      std::vector<double> & vfrac
      ){


    //if(m_is_verbose)           
    std::cout<<"RecoMCTracks()\t||\t Begininning recob::Track Reco-MC suite on: "<<tracks.size()<<" tracks."<<std::endl;

    int i_trk = 0;

    for(size_t k =0; k< tracks.size();++k){
      const art::Ptr<recob::Track> track = tracks[k];
      m_sim_track_matched[i_trk] = 0;

      if(trackToMCParticleMap.count(track)>0){

        const art::Ptr<simb::MCParticle> mcparticle = trackToMCParticleMap[track];
        std::cout<<"count2: "<<MCParticleToMCTruthMap.count(mcparticle)<<std::endl;
        const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruthMap[mcparticle];

        PandoraPFParticle* ppfp = PPFP_GetPPFPFromTrack(all_PPFPs, track);
        //        const art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;
        //    const art::Ptr<recob::PFParticle> pfp = //trackToPFParticleMap[track];

        std::vector<double> correctedstart(3);
        std::vector<double> correctedend(3);
        std::vector<double> raw_End  ={mcparticle->EndX(), mcparticle->EndY(), mcparticle->EndZ()};
        // std::cout<<"the raw end of this mcparticle is "<<raw_End[0]<<", "<<raw_End[1]<<", "<<raw_End[2]<<std::endl;
        spacecharge_correction(mcparticle, correctedstart);
        spacecharge_correction(mcparticle, correctedend, raw_End);

        //std::cout<<"the corrected end of this mcparticle is "<<correctedend[0]<<", "<<correctedend[1]<<", "<<correctedend[2]<<std::endl;


        m_sim_track_matched[i_trk] = 1;
        m_sim_track_energy[i_trk] = mcparticle->E();
        m_sim_track_mass[i_trk] = mcparticle->Mass();
        m_sim_track_kinetic_energy[i_trk] = m_sim_track_energy[i_trk]-m_sim_track_mass[i_trk];
        m_sim_track_pdg[i_trk] = mcparticle->PdgCode();
        m_sim_track_process[i_trk] = mcparticle->Process();
        m_sim_track_startx[i_trk] = correctedstart[0];
        m_sim_track_starty[i_trk] = correctedstart[1];
        m_sim_track_startz[i_trk] = correctedstart[2];

        m_sim_track_endx[i_trk]= correctedend[0];
        m_sim_track_endy[i_trk]= correctedend[1];
        m_sim_track_endz[i_trk]= correctedend[2];

        m_sim_track_length[i_trk]= sqrt(pow( m_sim_track_endx[i_trk] -  m_sim_track_startx[i_trk], 2)+ pow( m_sim_track_endy[i_trk] -  m_sim_track_starty[i_trk], 2) + pow( m_sim_track_endz[i_trk] -  m_sim_track_startz[i_trk], 2));

        m_sim_track_px[i_trk]=  mcparticle->Px();
        m_sim_track_py[i_trk]=  mcparticle->Py();
        m_sim_track_pz[i_trk]=  mcparticle->Pz();


        m_sim_track_origin[i_trk] = mctruth->Origin();
        m_sim_track_trackID[i_trk] = mcparticle->TrackId();
        m_sim_track_overlay_fraction[i_trk] = vfrac[i_trk];

        m_sim_track_sliceId[i_trk] = ppfp->get_SliceID();//PFPToSliceIdMap[pfp];
        m_sim_track_nuscore[i_trk] = ppfp->get_NuScore();//sliceIdToNuScoreMap[ m_sim_track_sliceId[i_trk]] ;
        m_sim_track_isclearcosmic[i_trk] = ppfp->get_IsClearCosmic();//PFPToClearCosmicMap[pfp]; 


        if(mcparticle->Mother()>=(int)mcParticleVector.size()){
          m_sim_track_parent_pdg[i_trk] = -1;
        }else{
          m_sim_track_parent_pdg[i_trk] = mcParticleVector[mcparticle->Mother()]->PdgCode();
        }

      }
      i_trk++;
    }

    return;
  }
  //recoMCmatching but specifically for recob::showers
  void showerRecoMCmatching(
      std::vector<PandoraPFParticle> all_PPFPs,
      std::vector<art::Ptr<recob::Shower>>& showerVector,
      std::map<art::Ptr<recob::Shower>,art::Ptr<simb::MCParticle>>& showerToMCParticleMap,
      art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
      std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
      std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap){

    std::vector<double> vec_fraction_matched;
    //processes that are "showery"
    std::map<std::string,bool> map_is_shower_process = {{"compt",true},
      {"FastScintillation",true},
      {"eBrem",true},
      {"phot",true},
      {"eIoni",true},
      {"conv",true},
      {"annihil",true}};

    std::vector<int> spacers = Printer_header({"pfpID", " matched_#simb", " pdg","    E_plane0","    E_plane1","    E_plane2"," cosmic?"});
    //for each recob::track/shower in the event
    bool default_verbose = m_is_verbose;
    m_is_verbose = false;
    for(size_t i=0; i<showerVector.size();++i){
      auto shower = showerVector[i];

      //get the associated reco PFP
      //            if(m_is_verbose)std::cout<<"We have "<<showerToPFParticleMap.count(shower)<<" matches in map"<<std::endl;

      PandoraPFParticle* ppfp = PPFP_GetPPFPFromShower(all_PPFPs, shower);
      const art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;

      //putting in the PFP pdg code as a check

      //and get the hits associated to the reco PFP
      std::vector< art::Ptr<recob::Hit> > obj_hits_ptrs = ppfp->pPFPHits; //pfParticleToHitsMap[pfp];

      /**
       * 
       * Loop over hits associated with the reco PFP to find MCParticles which contribute energy to the reco shower
       *
       **/

      std::unordered_map<int,double> objide; //map between the MCParticle track ID and the backtracker energy

      //energy for an MCParticle that comprises the most energy when sum over associated hits in PFP
      //total energy of the reco PFP taken from the sum of the hits associated to an MCParticle
      double maxe=-1, tote=0;                

      std::vector<double> total_energy_on_plane = {0.0,0.0,0.0};
      //simb::MCParticle const * best_matched_mcparticle = NULL; //pointer for the particle match we will calculate
      art::Ptr<simb::MCParticle> best_matched_mcparticle; //pointer for the MCParticle match we will calculate

      //    std::vector<simb::MCParticle const *> particle_vec;
      //    std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

      std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
      std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

      int n_associated_mcparticle_hits = 0;
      int n_not_associated_hits = 0;

      //this is the vector that will store the associated MC paritcles, as well as a MAP to the amount of energy associated
      std::vector<art::Ptr<simb::MCParticle>> asso_mcparticles_vec;
      std::map<art::Ptr<simb::MCParticle>, std::vector<double>> map_asso_mcparticles_energy;

      bool found_a_match = false;

      //std::cout<<"RecoMC()\t||\t On shower: "<<i<<" with pfp "<< pfp->Self() <<"and slice id "<<PFPToSliceIdMap[pfp]<<". This shower has "<<obj_hits_ptrs.size()<<" hits associated with it"<<std::endl;

      //loop only over hits associated to this reco PFP
      for(size_t i_h=0; i_h < obj_hits_ptrs.size(); ++i_h){

        int which_plane = (int)obj_hits_ptrs[i_h]->View();

        particle_vec.clear(); match_vec.clear(); //only store per hit

        //for the hit, fill the backtracker info 
        mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(), particle_vec, match_vec);
        // std::cout<<"for hit "<< i_h <<" particle_vec.size() = "<< particle_vec.size()<< " and match_vec.size() = "<< match_vec.size()<<std::endl; 


        //mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
        //the .key() gives us the index in the original collection
        //std::cout<<"REC: hit "<<i_h<<" has "<<particle_vec.size()<<" MCparticles assocaied: "<<std::endl;

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


      } // end loop over hits

      double fraction_num_hits_overlay = (double)n_not_associated_hits/(double)obj_hits_ptrs.size();

      if(m_is_verbose)std::cout << "recoMC()\t||\t On Object "<<i<<". The number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       
      if(m_is_verbose) std::cout<<"recoMC()\t||\t the fraction of hits from overlay is is "<<fraction_num_hits_overlay<<" ("<<n_not_associated_hits<<"/"<<obj_hits_ptrs.size()<<")"<<std::endl;


      if(n_associated_mcparticle_hits == 0){
        //This will only occur if the whole recob::PFParticle is PURELY associated with an overlay shower
        found_a_match =false;
        if(!found_a_match){
        }
        //Here we will fill every sim_shower_XXX variable with -999 or something like that 

        m_sim_shower_matched[i] = 0;
        m_sim_shower_energy[i] = -999;
        m_sim_shower_mass[i] = -999;
        m_sim_shower_kinetic_energy[i] = -999;
        m_sim_shower_pdg[i] = -999;
        m_sim_shower_trackID[i] = -999;
        m_sim_shower_process[i] = "overlay";
        m_sim_shower_end_process[i] = "overlay";
        m_sim_shower_parent_pdg[i] = -999;
        m_sim_shower_parent_trackID[i] = -999;
        m_sim_shower_vertex_x[i] = -9999;
        m_sim_shower_vertex_y[i] = -9999;
        m_sim_shower_vertex_z[i] = -9999;

        m_sim_shower_start_x[i] = -9999;
        m_sim_shower_start_y[i] = -9999;
        m_sim_shower_start_z[i] = -9999;
        m_sim_shower_px[i] = -9999;
        m_sim_shower_py[i] = -9999;
        m_sim_shower_pz[i] = -9999;

        m_sim_shower_is_true_shower[i] = -999;
        m_sim_shower_best_matched_plane[i] = -999;
        m_sim_shower_matched_energy_fraction_plane0[i] = -999;
        m_sim_shower_matched_energy_fraction_plane1[i] = -999;
        m_sim_shower_matched_energy_fraction_plane2[i] = -999;

        m_sim_shower_overlay_fraction[i] = fraction_num_hits_overlay;
        m_sim_shower_sliceId[i] = -999;
        m_sim_shower_nuscore[i] = -999;
        m_sim_shower_isclearcosmic[i] = -999;
        m_sim_shower_is_nuslice[i] = -999;


        continue;
      }


      /*  ********** if shower has been matched to MCParticle ************************* */

      /*
       *
       * Loop over each MCParticle associated to the reco shower to find the source particle
       *
       */

      std::map<int, art::Ptr<simb::MCParticle>> mother_MCP_map; //map between MCP track id and the source MCP

      std::vector<art::Ptr<simb::MCParticle>> marks_mother_vector;  // a vector of mother MCP
      std::map<art::Ptr<simb::MCParticle>, std::vector<double>> marks_mother_energy_fraction_map; // map of mother MCP and its energy on 3 planes

      int this_mcp_id = -1; //the track id for the current MCP in parent tree
      int last_mcp_id = -1; //the track id for the previous MCP in parent tree
      int i_mcp = 0;

      int num_bt_mothers =0;  // number of associated MCP that has mothers

      //m_is_verbose = false;
      //for each MCP that's associated to the reco shower
      for(auto mcp:asso_mcparticles_vec){

        if(m_is_verbose) std::cout<<"-----------------------------Start L1 Loop --------------------------------------------------"<<std::endl;
        if(m_is_verbose) std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<") Start by Looking at an MCP with pdg code "<<mcp->PdgCode()<<" and status code "<<mcp->StatusCode()<<" TrackID: "<<mcp->TrackId()<<std::endl;
        if(m_is_verbose) std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<") This MCP gave "<<   map_asso_mcparticles_energy[mcp][0] <<" | "<<map_asso_mcparticles_energy[mcp][1]<<" | "<<map_asso_mcparticles_energy[mcp][2]<<" energy to the recob::Object on each plane"<<std::endl;
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
            if(m_is_verbose)   std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<")  null pointer at id "<<this_mcp_id<<std::endl;
            this_mcp_id = last_mcp_id; //if invalid, move back a level to the previous MCP in parent tree and break the loop
            break;
          }

          //If primary particle will have process "primary"
          if(m_is_verbose)    std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<")  going up the tree at an MCP with track id  "<<this_mcp_id<<", pdg code "<<this_mcp->PdgCode()<<", and status code "<<this_mcp->StatusCode()<<" and Mother: "<<this_mcp->Mother()<<" Process: "<<this_mcp->Process()<<" EndProcess: "<<this_mcp->EndProcess()<<std::endl;

          //if it is a valid particle, iterate forward to the mother
          last_mcp_id = this_mcp_id;
          this_mcp_id =  this_mcp->Mother();

          //Check to see if this MCP was created in a "showery" process
          if(map_is_shower_process.count(this_mcp->Process()) > 0){
            //if it was, keep going, 

          }else if(this_mcp->Process()=="primary"){
            //if its primary, great! Note it.
            if(m_is_verbose)  std::cout<<"L1: Backtracked to primary! breaking"<<std::endl;
            this_mcp_id = last_mcp_id; //if invalid, move back a level to the previous MCP in parent tree and break the loop
            break;
          }else{
            if(m_is_verbose) std::cout<<"L1: Backtracked to a particle created in "<<this_mcp->EndProcess()<<"! breaking"<<std::endl;
            this_mcp_id = last_mcp_id; //if invalid, move back a level to the previous MCP in parent tree and break the loop
            break;
          }
        }

        //if the MCP at the top of the interaction chain has a valid track id store this in the mother map
        if (this_mcp_id >= 0){

          //Guanqun: this line here doesn't really cosider other break cases than finding primary particle
          if(m_is_verbose)   std::cout<<"L1: ("<<i<<" <-> "<<i_mcp<<") Storing the mother mother particle with track id "<<this_mcp_id<<" and pdg code "<<MCParticleToTrackIdMap[this_mcp_id]->PdgCode()<<" and status code "<<MCParticleToTrackIdMap[this_mcp_id]->StatusCode()<<std::endl;

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
          if(m_is_verbose)  std::cout<<"L1: error, the mother mother id was "<<this_mcp_id <<std::endl;

        }

        if(m_is_verbose)  std::cout<<"-----------------------------End L1 Loop --------------------------------------------------"<<std::endl;
        i_mcp++;
      }//for each MCParticle that's associated to a the recob::Shower

      //m_is_verbose = true;
      //there should be at least 1 mother MCP
      if(m_is_verbose)           std::cout<<"recoMC()\t||\t the number of source mother particles is "<<mother_MCP_map.size()<<" of which : "<<marks_mother_vector.size()<<" are unique!"<<std::endl;

      if(m_is_verbose)       std::cout<<"---------------------------- L2-------------------------------"<<std::endl;

      double best_mother_index = 0;
      double best_mother_energy = -9999;
      //            int best_mother_plane = -99;

      for(size_t p=0; p< marks_mother_vector.size(); p++){
        art::Ptr<simb::MCParticle> mother = marks_mother_vector[p];
        std::vector<double> mother_energy_recod = marks_mother_energy_fraction_map[mother];
        if(m_is_verbose)    std::cout<<"L2: Mother candidate "<<p<<" TrackID "<<mother->TrackId()<<" Process: "<<mother->Process()<<" EndProcess: "<<mother->EndProcess()<<std::endl;
        if(m_is_verbose)   std::cout<<"L2: Mother candidate "<<p<<" Energy "<<mother->E()<<" Reco'd Energy: "<<mother_energy_recod[0]<<" | "<<mother_energy_recod[1]<<" | "<<mother_energy_recod[2]<<" Fraction: ("<<mother_energy_recod[0]/(1000*mother->E())*100.0<<"% , "<<mother_energy_recod[1]/(1000*mother->E())*100.0<<"% , "<<mother_energy_recod[2]/(1000*mother->E())*100.0<<"% )"<<std::endl;

        if( mother_energy_recod[0] > best_mother_energy){
          best_mother_index = p;
          best_mother_energy = mother_energy_recod[0];
          //                    best_mother_plane = 0;
        }

        if( mother_energy_recod[1] > best_mother_energy){
          best_mother_index = p;
          best_mother_energy = mother_energy_recod[1];
          //                    best_mother_plane = 1;
        }

        if( mother_energy_recod[2] > best_mother_energy){
          best_mother_index = p;
          best_mother_energy = mother_energy_recod[2];
          //                    best_mother_plane = 2;
        }

      }



      // now have found the best mother of the shower
      if(m_is_verbose) std::cout<<"---------------------------- L2-------------------------------"<<std::endl;
      const art::Ptr<simb::MCParticle> match = marks_mother_vector[best_mother_index];

      std::vector<double> corrected_vertex(3), corrected_start(3);
      spacecharge_correction(match, corrected_vertex);


      if(match->PdgCode()==22){ // if it's a gamma
        std::vector<double> tmp  ={match->EndX(), match->EndY(), match->EndZ()};
        spacecharge_correction(match, corrected_start, tmp );
        m_sim_shower_is_true_shower[i] = 1;
      }else if(abs(match->PdgCode())==11){  // if it's e+/e-
        spacecharge_correction(match, corrected_start);
        m_sim_shower_is_true_shower[i] = 1;
      }else{
        corrected_start = {-999,-999,-999};
        m_sim_shower_is_true_shower[i] = 0;
      }

      art::Ptr<simb::MCParticle> match_mother = MCParticleToTrackIdMap[match->Mother()];

      if (match_mother.isNull()){
        m_sim_shower_parent_pdg[i] = -1;
        m_sim_shower_parent_trackID[i] = -1;

      }else{
        m_sim_shower_parent_pdg[i] = match_mother->PdgCode();
        m_sim_shower_parent_trackID[i] = match_mother->TrackId();
      }



      m_sim_shower_matched[i] = 1;
      m_sim_shower_energy[i] = match->E();
      m_sim_shower_mass[i] = match->Mass();
      m_sim_shower_kinetic_energy[i] = match->E()-match->Mass();
      m_sim_shower_pdg[i] = match->PdgCode();
      m_sim_shower_trackID[i] = match->TrackId();
      m_sim_shower_process[i] = match->Process();
      m_sim_shower_end_process[i] = match->EndProcess();
      m_sim_shower_vertex_x[i] = corrected_vertex[0];
      m_sim_shower_vertex_y[i] = corrected_vertex[1];
      m_sim_shower_vertex_z[i] =corrected_vertex[2];

      m_sim_shower_start_x[i] = corrected_start[0];
      m_sim_shower_start_y[i] = corrected_start[1];
      m_sim_shower_start_z[i] =corrected_start[2];

      m_sim_shower_px[i] = match->Px();
      m_sim_shower_py[i] = match->Py();
      m_sim_shower_pz[i] = match->Pz();

      // should've use 'best_mother_plane' here
      m_sim_shower_best_matched_plane[i] = best_mother_index;
      m_sim_shower_matched_energy_fraction_plane0[i] = marks_mother_energy_fraction_map[marks_mother_vector[best_mother_index]][0]/total_energy_on_plane[0];
      m_sim_shower_matched_energy_fraction_plane1[i] = marks_mother_energy_fraction_map[marks_mother_vector[best_mother_index]][1]/total_energy_on_plane[1];
      m_sim_shower_matched_energy_fraction_plane2[i] = marks_mother_energy_fraction_map[marks_mother_vector[best_mother_index]][2]/total_energy_on_plane[2];

      m_sim_shower_overlay_fraction[i] = fraction_num_hits_overlay;

      mcParticleVector.push_back(match);
      showerToMCParticleMap[shower] = mcParticleVector.back();

      m_sim_shower_sliceId[i] = ppfp->get_SliceID();//PFPToSliceIdMap[pfp];
      m_sim_shower_nuscore[i] = ppfp->get_NuScore();//sliceIdToNuScoreMap[ m_sim_shower_sliceId[i]] ;
      m_sim_shower_isclearcosmic[i] = ppfp->get_IsClearCosmic();//PFPToClearCosmicMap[pfp];
      m_sim_shower_is_nuslice[i] = ppfp->get_IsNuSlice();//PFPToNuSliceMap[pfp];

      //            if(marks_mother_vector.size()!=0){
      //                //if(m_is_verbose)  std::cout<<"recoMC()\t||\t The `BEST` mother is a "<<marks_mother_vector[best_mother_index]->PdgCode()<<" at "<<best_mother_index<<" on plane: "<<best_mother_plane<<std::endl;
      //    //std::vector<int> spacers = Printer_header({"pfpID", " matched_#simb", " pdg","    E_plane0","    E_plane1","    E_plane2"," cosmic?"});
      //                std::cout<<"recoMC()\t||\t The `BEST` mother is a "<<marks_mother_vector[best_mother_index]->PdgCode()<<" at "<<best_mother_index<<" on plane: "<<best_mother_plane<<std::endl;
      //                for(int l=0; l<3; l++){
      //                    std::cout<<"recoMC()\t||\t It represents "<<marks_mother_energy_fraction_map[marks_mother_vector[best_mother_index]][l]/total_energy_on_plane[l]*100.0<<"% of the energy on plane: "<<l<<" which is "<<total_energy_on_plane[l] <<std::endl;
      //                }
      //            }
      //
      //
      //            if (m_sim_shower_isclearcosmic[i]== false){
      //                std::cout<<"sim shower is matched to non-clear cosmic PFP "<<pfp->Self()<<std::endl;
      //            }
      Printer_content({std::to_string(pfp->Self()),
          std::to_string(best_mother_index),
          std::to_string(marks_mother_vector[best_mother_index]->PdgCode()),
          std::to_string(total_energy_on_plane[0]),
          std::to_string(total_energy_on_plane[1]),
          std::to_string(total_energy_on_plane[2]),
          std::to_string(ppfp->get_IsClearCosmic())
          },spacers);

      if (m_is_verbose)  std::cout<<"looking at pfp "<< pfp->Self()<<" with is matched to true particle with pdg  m_sim_shower_pdg[i]= "<<  m_sim_shower_pdg[i]<< ". is_nuslice = "<< m_sim_shower_is_nuslice[i]<<" in slice "<< m_sim_shower_sliceId[i]<<". The matched energy for this shower from mark's mother particle with pdg "<<marks_mother_vector[best_mother_index]->PdgCode()<< " is "<<m_sim_shower_matched_energy_fraction_plane0[i]<<"/"<<m_sim_shower_matched_energy_fraction_plane1[i]<<"/" <<m_sim_shower_matched_energy_fraction_plane2[i]<<std::endl;

    }//end vector loop.
    m_is_verbose = default_verbose;
  }//end showerRecoMCmatching



  /* @brief: a simpler MCmatching function for track and shower
   * @argument to be filled in function body:
   *     objectToMCParticleMap: map of object (track, shower) to its best-matching MCParticle
   *     mcParticleVector: a vector of best-matching MCParticle corresponding to objectVector
   * @return: a vector of fraction number, which is the fraction of unassociated hits in all reco hits of PFParticle
   */  
  //Typenamed for recob::Track and recob::Shower
    std::vector<double> trackRecoMCmatching(std::vector<art::Ptr<recob::Track>>& objectVector,
        std::map<art::Ptr<recob::Track>,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
        std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
        art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
        std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector){

      std::vector<double> trk_overlay_vec;
      std::vector<double> vec_fraction_matched;
      bool reco_verbose = false;
      //for each recob::track/shower in the event
      for(size_t i=0; i<objectVector.size();++i){
        auto object = objectVector[i];

        //get the associated reco PFP
        const art::Ptr<recob::PFParticle> pfp = objectToPFParticleMap[object];

        // std::cout<<"recoMCmatching()\t||\t looking for a track match to pfp"<< pfp->Self()<<std::endl;


        int pdg = pfp->PdgCode();
        //and get the hits associated to the reco PFP
        std::vector< art::Ptr<recob::Hit> > obj_hits_ptrs = pfParticleToHitsMap[pfp];

        std::unordered_map<int,double> objide; //map between the MCParticle track ID and the backtracker energy

        //energy for an MCParticle that comprises the most energy when sum over associated hits in PFP
        //total energy of the reco PFP taken from the sum of the hits associated to an MCParticle
        double maxe=-1, tote=0;                

        //simb::MCParticle const * best_matched_mcparticle = NULL; //pointer for the particle match we will calculate
        art::Ptr<simb::MCParticle> best_matched_mcparticle; //pointer for the MCParticle match we will calculate

        //    std::vector<simb::MCParticle const *> particle_vec;
        //    std::vector<anab::BackTrackerHitMatchingData const *> match_vec;

        std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
        std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing

        bool found_a_match = false;
        int n_associated_mcparticle_hits = 0;
        int n_not_associated_hits = 0;

        //    std::cout<<"REC: This object with pfp "<< pfp->Self() <<" in slice "<<PFPToSliceIdMap[pfp] <<" has "<<obj_hits_ptrs.size()<<" hits associated with it"<<std::endl;

        //loop only over hits associated to this reco PFP
        for(size_t i_h=0; i_h < obj_hits_ptrs.size(); ++i_h){

          particle_vec.clear(); match_vec.clear(); //only store per hit

          //for the hit, fill the backtracker info

          mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(), particle_vec, match_vec);
          //    std::cout<<"for hit "<< i_h <<" particle_vec.size() = "<< particle_vec.size()<< " and match_vec.size() = "<< match_vec.size()<<std::endl; 

          //mcparticles_per_hit.get(obj_hits_ptrs[i_h].key(),particle_vec,match_vec);
          //the .key() gives us the index in the original collection
          //std::cout<<"REC: hit "<<i_h<<" has "<<particle_vec.size()<<" MCparticles assocaied: "<<std::endl;

          //if there is an MCParticle associated to this hit

          //if there is an MCParticle associated to this hit
          if(particle_vec.size()>0) n_associated_mcparticle_hits++;

          if(particle_vec.size()==0) n_not_associated_hits++;


          //loop over MCparticles finding which is the MCparticle with most "energy" matched correctly
          //for each MCParticle associated with this hit
          for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
            //add the energy of the back tracked hit for this MCParticle to the track id for the MCParticle in the map
            objide[ particle_vec[i_p]->TrackId()] += match_vec[i_p]->energy; //store energy per track id

            //add the energy of the back tracked hit to the total energy for the PFP
            tote += match_vec[i_p]->energy; //calculate total energy deposited

            //want the MCParticle with the max total energy summed from the back tracker hit energy from hits in PFP
            if( objide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
              maxe = objide[ particle_vec[i_p]->TrackId() ];
              best_matched_mcparticle = particle_vec[i_p];
              found_a_match = true;
            }
          }//end loop over particles per hit
        }


        double fraction_num_hits_overlay = (double)n_not_associated_hits/(double)obj_hits_ptrs.size();

        trk_overlay_vec.push_back(fraction_num_hits_overlay);
        if(n_associated_mcparticle_hits == 0){
          //This will only occur if the whole recob::PFParticle is associated with an overlay object
          //std::cout<<fraction_num_hits_overlay<<std::endl;
        }//for each recob::track/shower in the event

        //std::cout << "recoMC()\t||\t the number of MCParticles associated with this PFP is "<<objide.size()<<std::endl;       

        if(found_a_match){
          mcParticleVector.push_back(best_matched_mcparticle);
          objectToMCParticleMap[object] = mcParticleVector.back();
        }else{
          // mcParticleVector.push_back(0);
        }
        vec_fraction_matched.push_back(maxe/tote);
        // if(m_is_verbose){
        //     std::cout << "recoMC()\t||\t the fracrion matched is "<<maxe/tote<<std::endl;
        // }


        if(!found_a_match){
          if(reco_verbose) std::cout << "recoMC()\t||\t NO MATCH NO MATCH (from my loop) for PFP with pdg  "<<pdg<<std::endl;
          if(reco_verbose)std::cout<<" count "<<objectToMCParticleMap.count(object)<<std::endl;
        }else{
          //if(reco_verbose)  
          std::cout << "recoMC()\t||\t Final Match (from my loop) for PFP with pdg "<<pdg<<" is " << best_matched_mcparticle->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
            << " pdg=" << best_matched_mcparticle->PdgCode()
            << " trkid=" << best_matched_mcparticle->TrackId()
            << " ke=" << best_matched_mcparticle->E()-best_matched_mcparticle->Mass()<< "\n";
        }

      }//end vector loop.
      //return vec_fraction_matched;
      return trk_overlay_vec;
    }


  int    photoNuclearTesting(std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector){


    for(auto &mcp: mcParticleVector){
      int pdg = mcp->PdgCode();
      std::string end_process  = mcp->EndProcess();
      int status =  mcp->StatusCode()         ; 


      if(pdg==22){
        std::cout<<"PHOTO: "<<status<<" "<<end_process<<std::endl;
      }
    }


    return 0;

  }

}//namespace end
