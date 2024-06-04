#include <numeric>

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TPrincipal.h"

#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/analyze_PandoraReco.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/fiducial_volume.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/TruncMean.h"
#include "sbncode/SinglePhotonAnalysis/Libraries/Processors.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_math.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_gadget.h"


namespace single_photon
{



  //Analyze Tracks
  void AnalyzeTracks(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Track>>& tracks,
      std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfParticleToSpacePointsMap, 
      std::map<int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap,
      std::map<int, double> &sliceIdToNuScoreMap,
      var_all& vars,
      para_all& paras){


    if(g_is_verbose) std::cout<<"AnalyzeTracks()\t||\t Starting recob::Track analysis"<<std::endl;;

    int i_trk=0;


    //const double adc2eU(5.1e-3);
    //const double adc2eV(5.2e-3);
    // const double adc2eW(5.4e-3);


    //    const double tau(theDetector->ElectronLifetime());


    //loop over each recob::Track
    for (TrackVector::const_iterator iter = tracks.begin(), iterEnd = tracks.end(); iter != iterEnd; ++iter)
    {

      const art::Ptr<recob::Track> track = *iter;
      //            const art::Ptr<recob::PFParticle> pfp = trackToNuPFParticleMap[track];
      PandoraPFParticle* ppfp = PPFP_GetPPFPFromTrack(all_PPFPs, track);
      const art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;

      const std::vector< art::Ptr<recob::SpacePoint> > trk_spacepoints = pfParticleToSpacePointsMap[pfp];
      const std::vector<art::Ptr<recob::Hit>> trk_hits  = ppfp->pPFPHits;//pfParticleToHitsMap[pfp];

      int m_trkid = track->ID();
      double tem_length = track->Length();
      auto tem_trk_dir = track->Direction(); // type of vars.m_trk_dir: a std::pair of two 3D vector
      //first: direction of first point, second: direction of the end of track

      if(g_is_verbose) std::cout<<"AnalyzeTracks()\t||\t On Track: "<<i_trk<<" with TrackID: "<<m_trkid<<" and length: "<<tem_length<<""<<std::endl;;

      vars.m_reco_track_calo_energy_plane0[i_trk] = CalcEShowerPlane(trk_hits, 0, paras);
      vars.m_reco_track_calo_energy_plane1[i_trk] = CalcEShowerPlane(trk_hits, 1, paras);
      vars.m_reco_track_calo_energy_plane2[i_trk] = CalcEShowerPlane(trk_hits, 2, paras);
      vars.m_reco_track_calo_energy_max[i_trk] = std::max( vars.m_reco_track_calo_energy_plane2[i_trk],  std::max(vars.m_reco_track_calo_energy_plane0[i_trk],vars.m_reco_track_calo_energy_plane1[i_trk]));

      vars.m_reco_track_num_spacepoints[i_trk] = (int)trk_spacepoints.size();


      vars.m_reco_track_startx[i_trk]= track->Start().X();
      vars.m_reco_track_starty[i_trk]= track->Start().Y();
      vars.m_reco_track_startz[i_trk]= track->Start().Z();

      vars.m_reco_track_length[i_trk] =tem_length;
      vars.m_reco_track_dirx[i_trk] = tem_trk_dir.first.X();
      vars.m_reco_track_diry[i_trk] = tem_trk_dir.first.Y();
      vars.m_reco_track_dirz[i_trk] = tem_trk_dir.first.Z();

      vars.m_reco_track_endx[i_trk] = track->End().X();
      vars.m_reco_track_endy[i_trk]= track->End().Y();
      vars.m_reco_track_endz[i_trk]= track->End().Z();

      std::vector<double> hend = {vars.m_reco_track_endx[i_trk],vars.m_reco_track_endy[i_trk],vars.m_reco_track_endz[i_trk]};
      std::vector<double> hstart = {vars.m_reco_track_startx[i_trk],vars.m_reco_track_starty[i_trk],vars.m_reco_track_startz[i_trk]};

      vars.m_reco_track_end_dist_to_active_TPC[i_trk] =   distToTPCActive(hend, paras);
      vars.m_reco_track_start_dist_to_active_TPC[i_trk] = distToTPCActive(hstart, paras);

      vars.m_reco_track_end_dist_to_CPA[i_trk] = distToCPA(hend, paras);
      vars.m_reco_track_start_dist_to_CPA[i_trk] = distToCPA(hstart, paras);

      vars.m_reco_track_end_in_SCB[i_trk] = distToSCB(vars.m_reco_track_end_dist_to_SCB[i_trk],hend, paras);
      vars.m_reco_track_start_in_SCB[i_trk] = distToSCB(vars.m_reco_track_start_dist_to_SCB[i_trk],hstart, paras);


      vars.m_reco_track_theta_yz[i_trk] = atan2(vars.m_reco_track_diry[i_trk],vars.m_reco_track_dirz[i_trk]);
      vars.m_reco_track_phi_yx[i_trk] = atan2(vars.m_reco_track_diry[i_trk],vars.m_reco_track_dirx[i_trk]);

      std::vector<double> tmp_trk_start = {vars.m_reco_track_startx[i_trk],vars.m_reco_track_starty[i_trk],vars.m_reco_track_startz[i_trk]};
      std::vector<double> tmp_trk_end = {vars.m_reco_track_endx[i_trk],vars.m_reco_track_endy[i_trk],vars.m_reco_track_endz[i_trk]};
      double max_dist_from_line = -9999999;

      vars.m_reco_track_spacepoint_chi[i_trk] = 0.0;
      //Principal componant analysis of SPACEPOINTS and not trajectory points. Will add a few things here in future.

      if(g_is_verbose) std::cout<<"AnalyzeTracks()\t||\t Beginining PCA analysis of track"<<std::endl;;

      TPrincipal* principal;
      principal = new TPrincipal(3,"ND");
      for(int x = 0; x < vars.m_reco_track_num_spacepoints[i_trk]; x++){
        // get the position of spacepoint in xyz
        std::vector<double> tmp_spacepoints = {trk_spacepoints[x]->XYZ()[0],trk_spacepoints[x]->XYZ()[1] , trk_spacepoints[x]->XYZ()[2]};
        principal->AddRow(&tmp_spacepoints[0]);

        // distance between track direction and spacepoint
        double dist = dist_line_point(tmp_trk_start,tmp_trk_end,tmp_spacepoints);
        if(dist> max_dist_from_line) max_dist_from_line = dist;
        vars.m_reco_track_spacepoint_chi[i_trk] += dist*dist;
      }
      vars.m_reco_track_spacepoint_max_dist[i_trk]= max_dist_from_line;

      principal->MakePrincipals();
      TVectorD * eigen = (TVectorD*) principal->GetEigenValues();

      vars.m_reco_track_spacepoint_principal0[i_trk]=(*eigen)(0);
      vars.m_reco_track_spacepoint_principal1[i_trk]=(*eigen)(1);
      vars.m_reco_track_spacepoint_principal2[i_trk]=(*eigen)(2);

      delete principal;

      if(g_is_verbose) std::cout<<"AnalyzeTracks()\t||\t Completed PCA analysis of track. Primary componant: "<<vars.m_reco_track_spacepoint_principal0.back()<<""<<std::endl;;

      //range based energy calculation assuming

      if(paras.s_run_pi0_filter){
        // assume this track is a proton track, get its energy
        vars.m_reco_track_proton_kinetic_energy[i_trk] = -9999;
      }else{
        //WARNING, extra input is needed for the proton track energy;
        //        vars.m_reco_track_proton_kinetic_energy[i_trk] = proton_length2energy_tgraph.Eval(tem_length)/1000.0; 
      }

      if(tem_length == 0.0) vars.m_reco_track_proton_kinetic_energy[i_trk]=0.0;

      // Dead Wire Approximity
      //      vars.m_reco_track_end_to_nearest_dead_wire_plane0[i_trk] = distanceToNearestDeadWire(0, vars.m_reco_track_endy[i_trk], vars.m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);
      //      vars.m_reco_track_end_to_nearest_dead_wire_plane1[i_trk] = distanceToNearestDeadWire(1, vars.m_reco_track_endy[i_trk], vars.m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);
      //      vars.m_reco_track_end_to_nearest_dead_wire_plane2[i_trk] = distanceToNearestDeadWire(2, vars.m_reco_track_endy[i_trk], vars.m_reco_track_endz[i_trk],geom,bad_channel_list_fixed_mcc9);

      vars.m_reco_track_sliceId[i_trk] = ppfp->get_SliceID();//PFPToSliceIdMap[pfp];
      // Guanqun: how do you make sure the sliceId is positive, not -1, as for cosmic tracks?
      // sliceIdToNuScoreMap seems to only have sliceID:nuScore pairs for these with actual nuScores.
      vars.m_reco_track_nuscore[i_trk] = sliceIdToNuScoreMap[ vars.m_reco_track_sliceId[i_trk]] ;
      vars.m_reco_track_isclearcosmic[i_trk] = ppfp->get_IsClearCosmic();//PFPToClearCosmicMap[pfp];

      //std::cout<<"checking track nuslice"<<std::endl;
      // std::cout<<"is nuslice for track with pfp "<<pfp->Self()<<" is: "<<PFPToNuSliceMap[pfp]<<std::endl;
      vars.m_reco_track_is_nuslice[i_trk] = ppfp->get_IsNuSlice();//PFPToNuSliceMap[pfp];

      vars.m_reco_track_num_daughters[i_trk] = pfp->NumDaughters();
      if(vars.m_reco_track_num_daughters[i_trk]>0){
        //currently just look at 1 daughter
        //                vars.m_reco_track_daughter_trackscore[i_trk] = PFPToTrackScoreMap[pfParticleMap[pfp->Daughters().front()]];
        int pfp_size = all_PPFPs.size();
        for(int index = 0; index < pfp_size; index++){
          if( (pfp->Daughters().front()) == all_PPFPs[index].pPFParticle->Self()){
            vars.m_reco_track_daughter_trackscore[i_trk] = all_PPFPs[index].get_TrackScore();
            break;
          }
        }

      }

      vars.m_reco_track_trackscore[i_trk] = ppfp->get_TrackScore();
      vars.m_reco_track_pfparticle_pdg[i_trk] = ppfp->get_PdgCode();

      //  vars.m_reco_track_trackscore[i_trk] = PFPToTrackScoreMap[pfp];
      //            if ( PFPToTrackScoreMap.find(pfp) != PFPToTrackScoreMap.end() ) {
      //                vars.m_reco_track_trackscore[i_trk] = PFPToTrackScoreMap[pfp];
      //                vars.m_reco_track_pfparticle_pdg[i_trk] = pfp->PdgCode();
      //            } else{
      //                vars.m_reco_track_trackscore[i_trk] = -999; 
      //                vars.m_reco_track_pfparticle_pdg[i_trk] = -999; 
      //            }

      //A loop over the trajectory points
      size_t const traj_size = track->CountValidPoints();
      vars.m_reco_track_num_trajpoints[i_trk] = (int)traj_size;

      for(unsigned int  p = 0; p < traj_size; ++p) {
        //recob::Track::TrajectoryPoint_t const & trajp = track->TrajectoryPoint(j);
        //recob::Track::Point_t const & pos = trajp.position;
        //recob::Track::Vector_t const & mom = trajp.momentum;

      }


      i_trk++;

    } // end of recob::Track loop

    //Lets sort and order the showers
    vars.m_reco_track_ordered_energy_index = sort_indexes(vars.m_reco_track_proton_kinetic_energy);
    vars.m_reco_track_ordered_displacement_index = sort_indexes(vars.m_reco_track_length);


    if(g_is_verbose) std::cout<<"AnalyzeTracks()\t||\t Finished."<<std::endl;;
  }



  void AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::vector<PandoraPFParticle> all_PPFPs,var_all& vars, para_all& paras){

    if(g_is_verbose) std::cout<<"CollectCalo(recob::Track)\t||\t Starting calo module for recob::Track"<<std::endl;;

    for(size_t i_trk = 0; i_trk<tracks.size(); ++i_trk){
      const art::Ptr<recob::Track>      track = tracks[i_trk];
      PandoraPFParticle* ppfp = PPFP_GetPPFPFromTrack(all_PPFPs, track);
      std::vector<art::Ptr<anab::Calorimetry>> Calos = ppfp->get_Calorimetries();

      if(Calos.size()!=3){
        std::cout<<"Singlephoton::Tracks\t||\tERROR!! ERROR!!! anab::Calorimetery vector is not of length 3!!! "<<Calos.size()<<". Skipping this track!!"<<std::endl;
        continue;
      }
      //if(trackToCaloMap[track].size()!=3){
      //  std::cout<<"Singlephoton::Tracks\t||\tERROR!! ERROR!!! anab::Calorimetery vector is not of length 3!!! "<<trackToCaloMap[track].size()<<". Skipping this track!!"<<std::endl;
      //  continue;
      //}

      const art::Ptr<anab::Calorimetry> calo_p0  = Calos[0];
      const art::Ptr<anab::Calorimetry> calo_p1  = Calos[1];
      const art::Ptr<anab::Calorimetry> calo_p2  = Calos[2];


      size_t calo_length_p0 = calo_p0->dEdx().size();
      size_t calo_length_p1 = calo_p1->dEdx().size();
      size_t calo_length_p2 = calo_p2->dEdx().size();

      TruncMean tm_p0;
      TruncMean tm_p1;
      TruncMean tm_p2;

      std::vector<double> trunc_dEdx_p0;
      std::vector<double> res_range_good_p0;
      std::vector<double> dEdx_good_p0;

      std::vector<double> trunc_dEdx_p1;
      std::vector<double> res_range_good_p1;
      std::vector<double> dEdx_good_p1;

      std::vector<double> trunc_dEdx_p2;
      std::vector<double> res_range_good_p2;
      std::vector<double> dEdx_good_p2;

      vars.m_reco_track_num_calo_hits_p0[i_trk] = (int)calo_length_p0;
      vars.m_reco_track_num_calo_hits_p1[i_trk] = (int)calo_length_p1;
      vars.m_reco_track_num_calo_hits_p2[i_trk] = (int)calo_length_p2;

      vars.m_reco_track_best_calo_plane[i_trk]=-1;

      // guanqun: vectors have been clear and resized, so probably not need to reset their values?
      vars.m_reco_track_good_calo_p0[i_trk] =  0;
      vars.m_reco_track_mean_dEdx_p0[i_trk] =  0.0;
      vars.m_reco_track_mean_dEdx_start_half_p0[i_trk] =  0.0;
      vars.m_reco_track_mean_dEdx_end_half_p0[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_p0[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk] =  0.0;
      vars.m_reco_track_trunc_PIDA_p0[i_trk] =  0.0;

      vars.m_reco_track_good_calo_p1[i_trk] =  0;
      vars.m_reco_track_mean_dEdx_p1[i_trk] =  0.0;
      vars.m_reco_track_mean_dEdx_start_half_p1[i_trk] =  0.0;
      vars.m_reco_track_mean_dEdx_end_half_p1[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_p1[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk] =  0.0;
      vars.m_reco_track_trunc_PIDA_p1[i_trk] =  0.0;

      vars.m_reco_track_good_calo_p2[i_trk] =  0;
      vars.m_reco_track_mean_dEdx_p2[i_trk] =  0.0;
      vars.m_reco_track_mean_dEdx_start_half_p2[i_trk] =  0.0;
      vars.m_reco_track_mean_dEdx_end_half_p2[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_p2[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk] =  0.0;
      vars.m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk] =  0.0;
      vars.m_reco_track_trunc_PIDA_p2[i_trk] =  0.0;



      //First off look over ALL points
      //--------------------------------- plane 0 ----------- Induction
      for (size_t k = 0; k < calo_length_p0; ++k) {
        double res_range =    calo_p0->ResidualRange()[k];  //ResidualRange() returns range from end of track
        double dEdx =         calo_p0->dEdx()[k];

        vars.m_reco_track_mean_dEdx_p0[i_trk] += dEdx;
        if(k <= calo_length_p0/2){ 
          vars.m_reco_track_mean_dEdx_start_half_p0[i_trk]+=dEdx;
        }else{
          vars.m_reco_track_mean_dEdx_end_half_p0[i_trk]+=dEdx;
        }

        bool is_sensible = dEdx < paras.s_track_calo_max_dEdx; 
        bool is_nan =dEdx != dEdx; // != has higher precedence than = 
        bool is_inf = std::isinf(dEdx);
        bool is_nonzero = dEdx> paras.s_track_calo_min_dEdx;

        if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length_p0-1){
          res_range_good_p0.push_back(res_range);
          dEdx_good_p0.push_back(dEdx);
        }

        //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
      }// End of first loop.

      vars.m_reco_track_good_calo_p0[i_trk] = 0;
      if(res_range_good_p0.size() >= paras.s_track_calo_min_dEdx_hits){
        vars.m_reco_track_good_calo_p0[i_trk] = res_range_good_p0.size();

        //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimuvars.m_dx in the res_range
        double tenth_track = std::max(res_range_good_p0.front(), res_range_good_p0.back())/paras.s_track_calo_trunc_fraction;
        double min_dx = 999;
        for(int j = res_range_good_p0.size()-1; j>1; j--){
          double dx = fabs(res_range_good_p0[j]-res_range_good_p0[j-1]);
          if(dx < min_dx) min_dx = dx;
        }
        double rad = std::max( min_dx*2, tenth_track); 

        //Calculate the residual range
        tm_p0.setRadius(rad);
        tm_p0.CalcTruncMeanProfile(res_range_good_p0,dEdx_good_p0, trunc_dEdx_p0);      

        double pida_sum_trunc=0.0;
        //Calculate the mean truncated mean dEdx
        for(size_t k=0; k< trunc_dEdx_p0.size(); k++){
          double dEdx = trunc_dEdx_p0[k];
          vars.m_reco_track_mean_trunc_dEdx_p0[i_trk] += dEdx;
          if(k <= trunc_dEdx_p0.size()/2){ 
            vars.m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk]+=dEdx;
          }else{
            vars.m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk]+=dEdx;
          }


          if(trunc_dEdx_p0[k] != trunc_dEdx_p0[k] || std::isinf(trunc_dEdx_p0[k]) || trunc_dEdx_p0[k]<0){
            std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx_p0[k]<<std::endl;
            std::cout<<"Vector Length : "<<trunc_dEdx_p0.size()<<std::endl;
            std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
            //for(int m=0; m<trunc_dEdx.size(); m++){
            //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
            //}
            std::cout<<"Using Radius: "<<rad<<std::endl;
            //exit(EXIT_FAILURE);
            vars.m_reco_track_good_calo_p0[i_trk] = 0; 
          }

          // dEdx/pow(residual_range, -0.42) is the constant A in residual range formula
          pida_sum_trunc += trunc_dEdx_p0[k]/(pow(res_range_good_p0[k],-0.42));
        }
        vars.m_reco_track_trunc_PIDA_p0[i_trk] = pida_sum_trunc;           
        vars.m_reco_track_resrange_p0[i_trk] = res_range_good_p0;
        vars.m_reco_track_trunc_dEdx_p0[i_trk] = trunc_dEdx_p0;
        vars.m_reco_track_dEdx_p0[i_trk] = dEdx_good_p0;

        //std::cout<<"the residual range at the start is "<<res_range_good[0]<<std::endl;
      }

      vars.m_reco_track_mean_dEdx_p0[i_trk]            *=1.0/((double)calo_length_p0);
      vars.m_reco_track_mean_dEdx_start_half_p0[i_trk] *=2.0/((double)calo_length_p0);
      vars.m_reco_track_mean_dEdx_end_half_p0[i_trk]   *=2.0/((double)calo_length_p0);
      vars.m_reco_track_mean_trunc_dEdx_p0[i_trk]            *=1.0/((double)trunc_dEdx_p0.size());
      vars.m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk] *=2.0/((double)trunc_dEdx_p0.size());
      vars.m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk]   *=2.0/((double)trunc_dEdx_p0.size());
      vars.m_reco_track_trunc_PIDA_p0[i_trk]  *=1.0/((double)trunc_dEdx_p0.size());

      //First off look over ALL points
      //--------------------------------- plane 1 ----------- Induction
      for (size_t k = 0; k < calo_length_p1; ++k) {
        double res_range =    calo_p1->ResidualRange()[k];
        double dEdx =         calo_p1->dEdx()[k];

        vars.m_reco_track_mean_dEdx_p1[i_trk] += dEdx;
        if(k <= calo_length_p1/2){ 
          vars.m_reco_track_mean_dEdx_start_half_p1[i_trk]+=dEdx;
        }else{
          vars.m_reco_track_mean_dEdx_end_half_p1[i_trk]+=dEdx;
        }

        bool is_sensible = dEdx < paras.s_track_calo_max_dEdx; 
        bool is_nan =dEdx != dEdx; 
        bool is_inf = std::isinf(dEdx);
        bool is_nonzero = dEdx> paras.s_track_calo_min_dEdx;

        if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length_p1-1){
          res_range_good_p1.push_back(res_range);
          dEdx_good_p1.push_back(dEdx);
        }

        //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
      }// End of first loop.

      vars.m_reco_track_good_calo_p1[i_trk] = 0;
      if(res_range_good_p1.size() >= paras.s_track_calo_min_dEdx_hits){
        vars.m_reco_track_good_calo_p1[i_trk] = res_range_good_p1.size();

        //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimuvars.m_dx in the res_range
        double tenth_track = std::max(res_range_good_p1.front(), res_range_good_p1.back())/paras.s_track_calo_trunc_fraction;
        double min_dx = 999;
        for(int j = res_range_good_p1.size()-1; j>1; j--){
          double dx = fabs(res_range_good_p1[j]-res_range_good_p1[j-1]);
          if(dx < min_dx) min_dx = dx;
        }
        double rad = std::max( min_dx*2, tenth_track); 

        //Calculate the residual range
        tm_p1.setRadius(rad);
        tm_p1.CalcTruncMeanProfile(res_range_good_p1,dEdx_good_p1, trunc_dEdx_p1);      

        double pida_sum_trunc=0.0;
        //Calculate the mean truncated mean dEdx
        for(size_t k=0; k< trunc_dEdx_p1.size(); k++){
          double dEdx = trunc_dEdx_p1[k];
          vars.m_reco_track_mean_trunc_dEdx_p1[i_trk] += dEdx;
          if(k <= trunc_dEdx_p1.size()/2){ 
            vars.m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk]+=dEdx;
          }else{
            vars.m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk]+=dEdx;
          }


          if(trunc_dEdx_p1[k] != trunc_dEdx_p1[k] || std::isinf(trunc_dEdx_p1[k]) || trunc_dEdx_p1[k]<0){
            std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx_p1[k]<<std::endl;
            std::cout<<"Vector Length : "<<trunc_dEdx_p1.size()<<std::endl;
            std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
            //for(int m=0; m<trunc_dEdx.size(); m++){
            //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
            //}
            std::cout<<"Using Radius: "<<rad<<std::endl;
            //exit(EXIT_FAILURE);
            vars.m_reco_track_good_calo_p1[i_trk] = 0; 
          }

          pida_sum_trunc += trunc_dEdx_p1[k]/(pow(res_range_good_p1[k],-0.42));
        }
        vars.m_reco_track_trunc_PIDA_p1[i_trk] = pida_sum_trunc;           
        vars.m_reco_track_resrange_p1[i_trk] = res_range_good_p1;
        vars.m_reco_track_trunc_dEdx_p1[i_trk] = trunc_dEdx_p1;
        vars.m_reco_track_dEdx_p1[i_trk] = dEdx_good_p1;

        //std::cout<<"the residual range at the start is "<<res_range_good[0]<<std::endl;
      }

      vars.m_reco_track_mean_dEdx_p1[i_trk]            *=1.0/((double)calo_length_p1);
      vars.m_reco_track_mean_dEdx_start_half_p1[i_trk] *=2.0/((double)calo_length_p1);
      vars.m_reco_track_mean_dEdx_end_half_p1[i_trk]   *=2.0/((double)calo_length_p1);
      vars.m_reco_track_mean_trunc_dEdx_p1[i_trk]            *=1.0/((double)trunc_dEdx_p1.size());
      vars.m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk] *=2.0/((double)trunc_dEdx_p1.size());
      vars.m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk]   *=2.0/((double)trunc_dEdx_p1.size());
      vars.m_reco_track_trunc_PIDA_p1[i_trk]  *=1.0/((double)trunc_dEdx_p1.size());

      //First off look over ALL points
      //--------------------------------- plane 2 ----------- Collection
      for (size_t k = 0; k < calo_length_p2; ++k) {
        double res_range =    calo_p2->ResidualRange()[k];
        double dEdx =         calo_p2->dEdx()[k];

        vars.m_reco_track_mean_dEdx_p2[i_trk] += dEdx;
        if(k <= calo_length_p2/2){ 
          vars.m_reco_track_mean_dEdx_start_half_p2[i_trk]+=dEdx;
        }else{
          vars.m_reco_track_mean_dEdx_end_half_p2[i_trk]+=dEdx;
        }

        bool is_sensible = dEdx < paras.s_track_calo_max_dEdx; 
        bool is_nan =dEdx != dEdx; 
        bool is_inf = std::isinf(dEdx);
        bool is_nonzero = dEdx> paras.s_track_calo_min_dEdx;

        if(is_sensible && !is_nan && !is_inf && is_nonzero && k != 0 && k != calo_length_p2-1){
          res_range_good_p2.push_back(res_range);
          dEdx_good_p2.push_back(dEdx);
        }

        //    std::cout<<"\t"<<k<<" "<<calo->dEdx()[k]<<" "<<calo->ResidualRange()[k]<<" "<< ""<<std::endl;;
      }// End of first loop.

      vars.m_reco_track_good_calo_p2[i_trk] = 0;
      if(res_range_good_p2.size() >= paras.s_track_calo_min_dEdx_hits){
        vars.m_reco_track_good_calo_p2[i_trk] = res_range_good_p2.size();

        //The radius we truncate over is going to be the max of either 1/frac of a track or 2x the minimuvars.m_dx in the res_range
        double tenth_track = std::max(res_range_good_p2.front(), res_range_good_p2.back())/paras.s_track_calo_trunc_fraction;
        double min_dx = 999;
        for(int j = res_range_good_p2.size()-1; j>1; j--){
          double dx = fabs(res_range_good_p2[j]-res_range_good_p2[j-1]);
          if(dx < min_dx) min_dx = dx;
        }
        double rad = std::max( min_dx*2, tenth_track); 

        //Calculate the residual range
        tm_p2.setRadius(rad);
        tm_p2.CalcTruncMeanProfile(res_range_good_p2,dEdx_good_p2, trunc_dEdx_p2);      

        double pida_sum_trunc=0.0;
        //Calculate the mean truncated mean dEdx
        for(size_t k=0; k< trunc_dEdx_p2.size(); k++){
          double dEdx = trunc_dEdx_p2[k];
          vars.m_reco_track_mean_trunc_dEdx_p2[i_trk] += dEdx;
          if(k <= trunc_dEdx_p2.size()/2){ 
            vars.m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk]+=dEdx;
          }else{
            vars.m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk]+=dEdx;
          }


          if(trunc_dEdx_p2[k] != trunc_dEdx_p2[k] || std::isinf(trunc_dEdx_p2[k]) || trunc_dEdx_p2[k]<0){
            std::cout<<"Truncated dedx is either inf or nan (or negative) @ "<<k<<" "<<trunc_dEdx_p2[k]<<std::endl;
            std::cout<<"Vector Length : "<<trunc_dEdx_p2.size()<<std::endl;
            std::cout<<"i\t range \t dedx \t trunc dedx"<<std::endl;
            //for(int m=0; m<trunc_dEdx.size(); m++){
            //    std::cout<<m<<"\t"<<c_resrange.at(m)<<"  "<<c_dEdx.at(m)<<"  "<<trunc_dEdx.at(m)<<std::endl;
            //}
            std::cout<<"Using Radius: "<<rad<<std::endl;
            //exit(EXIT_FAILURE);
            vars.m_reco_track_good_calo_p2[i_trk] = 0; 
          }

          pida_sum_trunc += trunc_dEdx_p2[k]/(pow(res_range_good_p2[k],-0.42));
        }
        vars.m_reco_track_trunc_PIDA_p2[i_trk] = pida_sum_trunc;           
        vars.m_reco_track_resrange_p2[i_trk] = res_range_good_p2;
        vars.m_reco_track_trunc_dEdx_p2[i_trk] = trunc_dEdx_p2;
        vars.m_reco_track_dEdx_p2[i_trk] = dEdx_good_p2;

        //std::cout<<"the residual range at the start is "<<res_range_good[0]<<std::endl;
      }

      vars.m_reco_track_mean_dEdx_p2[i_trk]            *=1.0/((double)calo_length_p2);
      vars.m_reco_track_mean_dEdx_start_half_p2[i_trk] *=2.0/((double)calo_length_p2);
      vars.m_reco_track_mean_dEdx_end_half_p2[i_trk]   *=2.0/((double)calo_length_p2);
      vars.m_reco_track_mean_trunc_dEdx_p2[i_trk]            *=1.0/((double)trunc_dEdx_p2.size());
      vars.m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk] *=2.0/((double)trunc_dEdx_p2.size());
      vars.m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk]   *=2.0/((double)trunc_dEdx_p2.size());
      vars.m_reco_track_trunc_PIDA_p2[i_trk]  *=1.0/((double)trunc_dEdx_p2.size());


      //************************ Now
      //lets pick one as a default?

      if(vars.m_reco_track_good_calo_p2[i_trk]!=0){
        vars.m_reco_track_best_calo_plane[i_trk] = 2;

        vars.m_reco_track_mean_dEdx_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_p2[i_trk];
        vars.m_reco_track_mean_dEdx_start_half_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_start_half_p2[i_trk];
        vars.m_reco_track_mean_dEdx_end_half_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_end_half_p2[i_trk];
        vars.m_reco_track_good_calo_best_plane[i_trk] = vars.m_reco_track_good_calo_p2[i_trk];
        vars.m_reco_track_trunc_dEdx_best_plane[i_trk] = vars.m_reco_track_trunc_dEdx_p2[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_p2[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_start_half_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_start_half_p2[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_end_half_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_end_half_p2[i_trk];
        vars.m_reco_track_trunc_PIDA_best_plane[i_trk] = vars.m_reco_track_trunc_PIDA_p2[i_trk];
        vars.m_reco_track_resrange_best_plane[i_trk] = vars.m_reco_track_resrange_p2[i_trk];
        vars.m_reco_track_dEdx_best_plane [i_trk] = vars.m_reco_track_dEdx_p2[i_trk];

      }else if(vars.m_reco_track_good_calo_p0[i_trk] > vars.m_reco_track_good_calo_p1[i_trk] ){
        vars.m_reco_track_best_calo_plane[i_trk] = 0;

        vars.m_reco_track_mean_dEdx_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_p0[i_trk];
        vars.m_reco_track_mean_dEdx_start_half_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_start_half_p0[i_trk];
        vars.m_reco_track_mean_dEdx_end_half_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_end_half_p0[i_trk];
        vars.m_reco_track_good_calo_best_plane[i_trk] = vars.m_reco_track_good_calo_p0[i_trk];
        vars.m_reco_track_trunc_dEdx_best_plane[i_trk] = vars.m_reco_track_trunc_dEdx_p0[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_p0[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_start_half_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_start_half_p0[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_end_half_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_end_half_p0[i_trk];
        vars.m_reco_track_trunc_PIDA_best_plane[i_trk] = vars.m_reco_track_trunc_PIDA_p0[i_trk];
        vars.m_reco_track_resrange_best_plane[i_trk] = vars.m_reco_track_resrange_p0[i_trk];
        vars.m_reco_track_dEdx_best_plane [i_trk] = vars.m_reco_track_dEdx_p0[i_trk];


      }else if(vars.m_reco_track_good_calo_p1[i_trk]!=0){
        vars.m_reco_track_best_calo_plane[i_trk] = 1;

        vars.m_reco_track_mean_dEdx_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_p1[i_trk];
        vars.m_reco_track_mean_dEdx_start_half_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_start_half_p1[i_trk];
        vars.m_reco_track_mean_dEdx_end_half_best_plane[i_trk] = vars.m_reco_track_mean_dEdx_end_half_p1[i_trk];
        vars.m_reco_track_good_calo_best_plane[i_trk] = vars.m_reco_track_good_calo_p1[i_trk];
        vars.m_reco_track_trunc_dEdx_best_plane[i_trk] = vars.m_reco_track_trunc_dEdx_p1[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_p1[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_start_half_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_start_half_p1[i_trk];
        vars.m_reco_track_mean_trunc_dEdx_end_half_best_plane[i_trk] = vars.m_reco_track_mean_trunc_dEdx_end_half_p1[i_trk];
        vars.m_reco_track_trunc_PIDA_best_plane[i_trk] = vars.m_reco_track_trunc_PIDA_p1[i_trk];
        vars.m_reco_track_resrange_best_plane[i_trk] = vars.m_reco_track_resrange_p1[i_trk];
        vars.m_reco_track_dEdx_best_plane [i_trk] = vars.m_reco_track_dEdx_p1[i_trk];



      }else{
        vars.m_reco_track_best_calo_plane[i_trk] = -1; 
      }


    }
  }



  void CollectPID(std::vector<art::Ptr<recob::Track>> & tracks,
      std::vector<PandoraPFParticle> all_PPFPs,var_all& vars){

    for(size_t i_trk=0; i_trk<tracks.size(); ++i_trk){
      art::Ptr<recob::Track> track = tracks[i_trk];

      PandoraPFParticle* ppfp = PPFP_GetPPFPFromTrack(all_PPFPs, track);
      art::Ptr<anab::ParticleID> pid = ppfp->get_ParticleID();//trackToPIDMap[track];
      if (!ppfp->get_HasPID()) {
        std::cout << "[analyze_Tracks] bad PID object" << std::endl;
        continue;
      }

      // For each PID object, create vector of PID scores for each algorithm
      // Loop over this and get scores for algorithm of choice
      // But first, prepare garbage values, just in case
      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = pid->ParticleIDAlgScores();
      double pidScore_BL_mu_plane0 = -999;
      double pidScore_BL_mu_plane1 = -999;
      double pidScore_BL_mu_plane2 = -999;
      double pidScore_BL_p_plane0 = -999;
      double pidScore_BL_p_plane1 = -999;
      double pidScore_BL_p_plane2 = -999;
      double pidScore_BL_mip_plane0 = -999;
      double pidScore_BL_mip_plane1 = -999;
      double pidScore_BL_mip_plane2 = -999;
      double pidScore_PIDA_plane0 = -999;
      double pidScore_PIDA_plane1 = -999;
      double pidScore_PIDA_plane2 = -999;
      double pidScore_chi2_mu_plane0 = -999;
      double pidScore_chi2_mu_plane1 = -999;
      double pidScore_chi2_mu_plane2 = -999;
      double pidScore_chi2_p_plane0 = -999;
      double pidScore_chi2_p_plane1 = -999;
      double pidScore_chi2_p_plane2 = -999;
      double pidScore_three_plane_proton = -999;

      //int planeid = 2;
      for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++) {

        anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
        //int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);
        int planeid = 0;//CHECK UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
        std::cout<<"WARNING: planeid is not set properly; so please check if you think this section matters"<<std::endl;
        /*
           if (planeid != 2){
           std::cout << "[ParticleIDValidation] Not using information for plane " 
           << planeid << " (using plane 2 calorimetry only)" << std::endl;
           continue;
           }
           */
        if (AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
          if (TMath::Abs(AlgScore.fAssumedPdg == 13)) {
            if (planeid == 0) pidScore_BL_mu_plane0 = AlgScore.fValue; //AlgScore.fValue: value produced by algorithm
            if (planeid == 1) pidScore_BL_mu_plane1 = AlgScore.fValue;
            if (planeid == 2) pidScore_BL_mu_plane2 = AlgScore.fValue;
          }
          if (TMath::Abs(AlgScore.fAssumedPdg == 2212)) {
            if (planeid == 0) pidScore_BL_p_plane0 = AlgScore.fValue;
            if (planeid == 1) pidScore_BL_p_plane1 = AlgScore.fValue;
            if (planeid == 2) pidScore_BL_p_plane2 = AlgScore.fValue;
          }
          if (TMath::Abs(AlgScore.fAssumedPdg == 0)) {
            if (planeid == 0) pidScore_BL_mip_plane0 = AlgScore.fValue;
            if (planeid == 1) pidScore_BL_mip_plane1 = AlgScore.fValue;
            if (planeid == 2) pidScore_BL_mip_plane2 = AlgScore.fValue;
          }
        }
        if (AlgScore.fAlgName == "Chi2" && anab::kVariableType(AlgScore.fVariableType) == anab::kGOF){
          if (TMath::Abs(AlgScore.fAssumedPdg == 13)) {
            if (planeid == 0) pidScore_chi2_mu_plane0 = AlgScore.fValue;
            if (planeid == 1) pidScore_chi2_mu_plane1 = AlgScore.fValue;
            if (planeid == 2) pidScore_chi2_mu_plane2 = AlgScore.fValue;
          }
          if (TMath::Abs(AlgScore.fAssumedPdg == 2212)) {
            if (planeid == 0) pidScore_chi2_p_plane0 = AlgScore.fValue;
            if (planeid == 1) pidScore_chi2_p_plane1 = AlgScore.fValue;
            if (planeid == 2) pidScore_chi2_p_plane2 = AlgScore.fValue;
          }
        }
        if (AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
          if (planeid == 0) pidScore_PIDA_plane0 = AlgScore.fValue;
          if (planeid == 1) pidScore_PIDA_plane1 = AlgScore.fValue;
          if (planeid == 2) pidScore_PIDA_plane2 = AlgScore.fValue;
        }
        //CHECK                    if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fPlaneMask == UBPID::uB_SinglePlaneGetBitset(2) ){}
        if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fPlaneMask == 0 ){
          //if (AlgScore.fAlgName == "ThreePlaneProtonPID" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood && TMath::Abs(AlgScore.fAssumedPdg) == 2212 && AlgScore.fPlaneMask == UBPID::uB_SinglePlaneGetBitset(2) && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){}
          pidScore_three_plane_proton = AlgScore.fValue;
        }
      }  // end looping over AlgScoresVec

      vars.m_reco_track_pid_bragg_likelihood_mu_plane0[i_trk] = pidScore_BL_mu_plane0;
      vars.m_reco_track_pid_bragg_likelihood_mu_plane1[i_trk] = pidScore_BL_mu_plane1;
      vars.m_reco_track_pid_bragg_likelihood_mu_plane2[i_trk] = pidScore_BL_mu_plane2;
      vars.m_reco_track_pid_bragg_likelihood_p_plane0[i_trk] = pidScore_BL_p_plane0;
      vars.m_reco_track_pid_bragg_likelihood_p_plane1[i_trk] = pidScore_BL_p_plane1;
      vars.m_reco_track_pid_bragg_likelihood_p_plane2[i_trk] = pidScore_BL_p_plane2;
      vars.m_reco_track_pid_bragg_likelihood_mip_plane0[i_trk] = pidScore_BL_mip_plane0;
      vars.m_reco_track_pid_bragg_likelihood_mip_plane1[i_trk] = pidScore_BL_mip_plane1;
      vars.m_reco_track_pid_bragg_likelihood_mip_plane2[i_trk] = pidScore_BL_mip_plane2;
      vars.m_reco_track_pid_chi2_mu_plane0[i_trk] = pidScore_chi2_mu_plane0;
      vars.m_reco_track_pid_chi2_mu_plane1[i_trk] = pidScore_chi2_mu_plane1;
      vars.m_reco_track_pid_chi2_mu_plane2[i_trk] = pidScore_chi2_mu_plane2;
      vars.m_reco_track_pid_chi2_p_plane0[i_trk] = pidScore_chi2_p_plane0;
      vars.m_reco_track_pid_chi2_p_plane1[i_trk] = pidScore_chi2_p_plane1;
      vars.m_reco_track_pid_chi2_p_plane2[i_trk] = pidScore_chi2_p_plane2;
      vars.m_reco_track_pid_pida_plane0[i_trk] = pidScore_PIDA_plane0;
      vars.m_reco_track_pid_pida_plane1[i_trk] = pidScore_PIDA_plane1;
      vars.m_reco_track_pid_pida_plane2[i_trk] = pidScore_PIDA_plane2;
      vars.m_reco_track_pid_three_plane_proton_pid[i_trk] = pidScore_three_plane_proton;

    }
    return;
  }



  //Analyze falshes
  void AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes, art::Handle<std::vector<sbn::crt::CRTHit>> crthit_h, double evt_timeGPS_nsec,  std::map<art::Ptr<recob::OpFlash>, std::vector< art::Ptr<sbn::crt::CRTHit>>> crtvetoToFlashMap,
      var_all& vars, para_all& paras){


    for(auto pair: crtvetoToFlashMap){
      std::cout<<"for flash at time "<< pair.first->Time()<<" has "<< pair.second.size() << " associated  CRT hits "<<std::endl;
      if(pair.second.size() > 0){
        for (auto hit: pair.second){
          std::cout<<"---- associated CRT hit at time "<<hit->ts0_ns/1000. <<" with PE "<<hit->peshit<<std::endl;
          vars.m_CRT_veto_hit_PE.push_back(hit->peshit);
        }

      }
      vars.m_CRT_veto_nhits =  pair.second.size();//save the number of associated CRT veto hits
    }


    if(g_is_verbose) std::cout<<"AnalyzeFlashes()\t||\t Beginning analysis of recob::OpFlash\n";

    size_t flash_size = flashes.size();
    for(size_t i = 0; i < flash_size; ++i) {

      art::Ptr<recob::OpFlash> const & flash = flashes[i];

      vars.m_reco_flash_total_pe[i]=(flash->TotalPE());
      vars.m_reco_flash_time[i]=(flash->Time());
      vars.m_reco_flash_time_width[i]=flash->TimeWidth();
      vars.m_reco_flash_abs_time[i]=flash->AbsTime();
      vars.m_reco_flash_frame[i]=flash->Frame();
      vars.m_reco_flash_ycenter[i]=flash->YCenter();
      vars.m_reco_flash_ywidth[i]=flash->YWidth();
      vars.m_reco_flash_zcenter[i]=flash->ZCenter();
      vars.m_reco_flash_zwidth[i]=flash->ZWidth();

      // paras.s_beamgate_flash_end/paras.s_beamgate_flash_start are read from pset
      if(vars.m_reco_flash_time[i] <= paras.s_beamgate_flash_end && vars.m_reco_flash_time[i] >= paras.s_beamgate_flash_start){
        vars.m_reco_num_flashes_in_beamgate++;
        vars.m_reco_flash_total_pe_in_beamgate[i]=(flash->TotalPE());
        vars.m_reco_flash_time_in_beamgate[i]=(flash->Time());
        vars.m_reco_flash_ycenter_in_beamgate[i] = flash->YCenter();
        vars.m_reco_flash_zcenter_in_beamgate[i] = flash->ZCenter();
      }

    }

    if(g_is_verbose) std::cout<<"AnalyzeFlashes()\t||\t Finished. There was "<<flash_size<<" flashes with: "<<vars.m_reco_num_flashes_in_beamgate<<" in the beamgate defined by: "<<paras.s_beamgate_flash_start<<" <-> "<<paras.s_beamgate_flash_end<<std::endl;

    //fill these values only for events that have CRT information - run3 G and later
    //code taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter_module.cc
    if(paras.s_runCRT){
      if (vars.m_reco_num_flashes_in_beamgate == 1){ //fill only if there's a flash in the beamgate

        int  _nCRThits_in_event = crthit_h->size();

        double _dt_abs   = 100000.0;
        //  double  _within_resolution = 0;
        double _beam_flash_time  =  vars.m_reco_flash_time_in_beamgate[0];  // Guanqun: why use index 0?

        // Loop over the CRT hits.
        for (int j = 0; j < _nCRThits_in_event; j++)
        {
          /*
             if (verbose)
             std::cout << "\t Time of the CRT Hit wrt the event timestamp = " << ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.) << " us." << std::endl;
             */
          double _crt_time_temp = ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + paras.s_DTOffset) / 1000.);

          // Fill the vector variables.
          vars.m_CRT_hits_time.push_back(_crt_time_temp);
          vars.m_CRT_hits_PE.push_back(crthit_h->at(j).peshit);
          vars.m_CRT_hits_x.push_back(crthit_h->at(j).x_pos);
          vars.m_CRT_hits_y.push_back(crthit_h->at(j).y_pos);
          vars.m_CRT_hits_z.push_back(crthit_h->at(j).z_pos);

          if (fabs(_beam_flash_time - _crt_time_temp) < _dt_abs)
          {
            _dt_abs = fabs(_beam_flash_time - _crt_time_temp);
            vars.m_CRT_dt = _beam_flash_time - _crt_time_temp;
            vars.m_CRT_min_hit_time = _crt_time_temp;
            // set 'within_resolution' to 'true' and break the loop if 'closest_crt_diff' is less than fResolution.
            if (_dt_abs < paras.s_Resolution)
            {
              //_within_resolution = 1;
              // Set the position information and the intensity of the CRT hit.
              vars.m_CRT_min_hit_PE = crthit_h->at(j).peshit;
              vars.m_CRT_min_hit_x = crthit_h->at(j).x_pos;
              vars.m_CRT_min_hit_y = crthit_h->at(j).y_pos;
              vars.m_CRT_min_hit_z = crthit_h->at(j).z_pos;


              // if (verbose)
              // {
              std::cout << "CRT hit PE = " << vars.m_CRT_min_hit_PE << " PEs." << std::endl;
              std::cout << "CRT hit x = " << vars.m_CRT_min_hit_x << " cm." << std::endl;
              std::cout << "CRT hit y = " << vars.m_CRT_min_hit_y << " cm." << std::endl;
              std::cout << "CRT hit z = " << vars.m_CRT_min_hit_z << " cm." << std::endl;
              // }
              break;
            }
          } // End of conditional for closest CRT hit time.
        } // End of loop over CRT hits.
      } //if there is 1 flash in beamgate
    }//if runCRT
  }//analyze flashes




  //Analyze Showers
  void AnalyzeShowers(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Shower>>& showers,  
      std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> >  & clusterToHitMap , 
      double triggeroffset,
      detinfo::DetectorPropertiesData const & theDetector,
      var_all& vars,
      para_all& paras){
    //        if(g_is_verbose) std::cout<<"AnalyzeShowers()\t||\t Begininning recob::Shower analysis suite"<<std::endl;;

    int i_shr = 0;

    std::vector<int> spacers = Printer_header({"Slice"," pfpID"," Start(x,    ","   y,      ","       z  )"," trackscore"," pdg"});
    for (ShowerVector::const_iterator iter = showers.begin(), iterEnd = showers.end(); iter != iterEnd; ++iter)
    {


      const art::Ptr<recob::Shower> shower = *iter;
      //            const art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap[shower];
      PandoraPFParticle* ppfp = PPFP_GetPPFPFromShower(all_PPFPs, shower);

      const art::Ptr<recob::PFParticle> pfp = ppfp->pPFParticle;

      art::Ptr<recob::Shower> shower3d;
      vars.m_reco_shower3d_exists[i_shr] = 0;
      shower3d = shower;

      const std::vector<art::Ptr<recob::Hit>> hits =  ppfp->pPFPHits;
      const std::vector<art::Ptr<recob::Cluster>> clusters = ppfp->pClusters;

      //int vars.m_shrid = shower->ID(); This is an used variable, always -999
      double tem_length = shower->Length();
      double tem_open_angle = shower->OpenAngle();

      TVector3 shr_start = shower->ShowerStart();
      TVector3 shr_dir = shower->Direction();

      TVector3 shr3d_start = shower3d->ShowerStart();
      TVector3 shr3d_dir = shower3d->Direction();

      //            if(g_is_verbose) std::cout<<"AnalyzeShowers()\t||\t On Shower: "<<i_shr<<" which has length: "<<tem_length<<""<<std::endl;;

      vars.m_reco_shower_startx[i_shr] = shr_start.X();
      vars.m_reco_shower_starty[i_shr] = shr_start.Y();
      vars.m_reco_shower_startz[i_shr] = shr_start.Z();


      std::vector<double> hstart = {vars.m_reco_shower_startx[i_shr],vars.m_reco_shower_starty[i_shr],vars.m_reco_shower_startz[i_shr]};
      vars.m_reco_shower_start_dist_to_active_TPC[i_shr] = distToTPCActive(hstart, paras);
      vars.m_reco_shower_start_dist_to_CPA[i_shr]        = distToCPA(hstart, paras);
      vars.m_reco_shower_start_in_SCB[i_shr] = distToSCB(vars.m_reco_shower_start_dist_to_SCB[i_shr],hstart, paras);

      vars.m_reco_shower_dirx[i_shr] = shr_dir.X();
      vars.m_reco_shower_diry[i_shr] = shr_dir.Y();
      vars.m_reco_shower_dirz[i_shr] = shr_dir.Z();
      vars.m_reco_shower_length[i_shr] = tem_length;
      vars.m_reco_shower_openingangle[i_shr] = tem_open_angle;

      vars.m_reco_shower3d_startx[i_shr] = shr3d_start.X();
      vars.m_reco_shower3d_starty[i_shr] = shr3d_start.Y();
      vars.m_reco_shower3d_startz[i_shr] = shr3d_start.Z();
      vars.m_reco_shower3d_dirx[i_shr] = shr3d_dir.X();
      vars.m_reco_shower3d_diry[i_shr] = shr3d_dir.Y();
      vars.m_reco_shower3d_dirz[i_shr] = shr3d_dir.Z();
      vars.m_reco_shower3d_length[i_shr] = shower3d->Length();
      vars.m_reco_shower3d_openingangle[i_shr] = shower3d->OpenAngle();


      vars.m_reco_shower_conversion_distance[i_shr] = sqrt( pow(shr_start.X()-vars.m_vertex_pos_x,2)+pow(shr_start.Y()-vars.m_vertex_pos_y,2)+ pow(shr_start.Z()-vars.m_vertex_pos_z,2)  );
      vars.m_reco_shower3d_conversion_distance[i_shr] = sqrt( pow(shr3d_start.X()-vars.m_vertex_pos_x,2)+pow(shr3d_start.Y()-vars.m_vertex_pos_y,2)+ pow(shr3d_start.Z()-vars.m_vertex_pos_z,2)  );

      //pandroa shower
      std::vector<double> shr_ts = {shr_start.X(), shr_start.Y(), shr_start.Z()};
      std::vector<double> shr_te = {shr_start.X()-shr_dir.X(),shr_start.Y()-shr_dir.Y(),shr_start.Z()-shr_dir.Z()};
      std::vector<double> shr_tv = {vars.m_vertex_pos_x,vars.m_vertex_pos_y,vars.m_vertex_pos_z};

      vars.m_reco_shower_impact_parameter[i_shr] = dist_line_point(shr_ts,shr_te,shr_tv );
      vars.m_reco_shower_implied_dirx[i_shr] = shr_start.X()-vars.m_vertex_pos_x;;
      vars.m_reco_shower_implied_diry[i_shr] = shr_start.Y()-vars.m_vertex_pos_y;
      vars.m_reco_shower_implied_dirz[i_shr] = shr_start.Z()-vars.m_vertex_pos_z;

      double norm = sqrt(pow(vars.m_reco_shower_implied_dirx[i_shr],2)+pow(vars.m_reco_shower_implied_diry[i_shr],2)+pow(vars.m_reco_shower_implied_dirz[i_shr],2));
      vars.m_reco_shower_implied_dirx[i_shr] = vars.m_reco_shower_implied_dirx[i_shr]/norm;
      vars.m_reco_shower_implied_diry[i_shr] = vars.m_reco_shower_implied_diry[i_shr]/norm;
      vars.m_reco_shower_implied_dirz[i_shr] = vars.m_reco_shower_implied_dirz[i_shr]/norm;

      //now 3D shower
      std::vector<double> shr3d_ts = {shr3d_start.X(), shr3d_start.Y(), shr3d_start.Z()};
      std::vector<double> shr3d_te = {shr3d_start.X()-shr3d_dir.X(),shr3d_start.Y()-shr3d_dir.Y(),shr3d_start.Z()-shr3d_dir.Z()};
      std::vector<double> shr3d_tv = {vars.m_vertex_pos_x,vars.m_vertex_pos_y,vars.m_vertex_pos_z};

      vars.m_reco_shower3d_impact_parameter[i_shr] = dist_line_point(shr3d_ts,shr3d_te,shr3d_tv );
      vars.m_reco_shower3d_implied_dirx[i_shr] = shr3d_start.X()-vars.m_vertex_pos_x;;
      vars.m_reco_shower3d_implied_diry[i_shr] = shr3d_start.Y()-vars.m_vertex_pos_y;
      vars.m_reco_shower3d_implied_dirz[i_shr] = shr3d_start.Z()-vars.m_vertex_pos_z;

      double shr3d_norm = sqrt(pow(vars.m_reco_shower3d_implied_dirx[i_shr],2)+pow(vars.m_reco_shower3d_implied_diry[i_shr],2)+pow(vars.m_reco_shower3d_implied_dirz[i_shr],2));
      vars.m_reco_shower3d_implied_dirx[i_shr] = vars.m_reco_shower3d_implied_dirx[i_shr]/shr3d_norm;
      vars.m_reco_shower3d_implied_diry[i_shr] = vars.m_reco_shower3d_implied_diry[i_shr]/shr3d_norm;
      vars.m_reco_shower3d_implied_dirz[i_shr] = vars.m_reco_shower3d_implied_dirz[i_shr]/shr3d_norm;


      vars.m_reco_shower_theta_yz[i_shr] = atan2(vars.m_reco_shower_diry[i_shr],vars.m_reco_shower_dirz[i_shr]);
      vars.m_reco_shower_phi_yx[i_shr] = atan2(vars.m_reco_shower_diry[i_shr],vars.m_reco_shower_dirx[i_shr]);

      vars.m_reco_shower3d_theta_yz[i_shr] = atan2(vars.m_reco_shower3d_diry[i_shr],vars.m_reco_shower3d_dirz[i_shr]);
      vars.m_reco_shower3d_phi_yx[i_shr] = atan2(vars.m_reco_shower3d_diry[i_shr],vars.m_reco_shower3d_dirx[i_shr]);


      //      vars.m_reco_shower_start_to_nearest_dead_wire_plane0[i_shr] = distanceToNearestDeadWire(0, vars.m_reco_shower_starty[i_shr], vars.m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
      //      vars.m_reco_shower_start_to_nearest_dead_wire_plane1[i_shr] = distanceToNearestDeadWire(1, vars.m_reco_shower_starty[i_shr], vars.m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
      //      vars.m_reco_shower_start_to_nearest_dead_wire_plane2[i_shr] = distanceToNearestDeadWire(2, vars.m_reco_shower_starty[i_shr], vars.m_reco_shower_startz[i_shr],geom, bad_channel_list_fixed_mcc9);
      std::vector<int> t_num(3,0);   // num of triangles on each plane
      std::vector<int> t_numhits(3,0);  // num of hits on each plane
      std::vector<double> t_area(3,0.0);

      //Right, this basically loops over all hits in all planes and for each plane forms the Delaunay triangilization of it and calculates the 2D area inscribed by the convex hull
      //            if(g_is_verbose) std::cout<<"AnalyzeShowers()\t||\t Starting Delaunay Triangleization"<<std::endl;;

      //auto start = std::chrono::high_resolution_clock::now();
      delaunay_hit_wrapper(hits, t_numhits, t_num, t_area, paras);

      //auto finish = std::chrono::high_resolution_clock::now();
      //auto microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
      //if(g_is_verbose) std::cout<<"AnalyzeShowers()\t||\t Finished Delaunay Triangleization. It took "<< microseconds.count() << "ms and found "<<t_num[0]+t_num[1]+t_num[2]<<" triangles"<<std::endl;;

      vars.m_reco_shower_delaunay_num_triangles_plane0[i_shr] = t_num[0];
      vars.m_reco_shower_delaunay_num_triangles_plane1[i_shr] = t_num[1];
      vars.m_reco_shower_delaunay_num_triangles_plane2[i_shr] = t_num[2];

      vars.m_reco_shower_delaunay_area_plane0[i_shr] = t_area[0];
      vars.m_reco_shower_delaunay_area_plane1[i_shr] = t_area[1];
      vars.m_reco_shower_delaunay_area_plane2[i_shr] = t_area[2];

      vars.m_reco_shower_num_hits_plane0[i_shr] = t_numhits[0];
      vars.m_reco_shower_num_hits_plane1[i_shr] = t_numhits[1];
      vars.m_reco_shower_num_hits_plane2[i_shr] = t_numhits[2];
      //-------------- Calorimetry 3D --------------------


      const std::vector< double > shr3d_energy = shower3d->Energy();
      const std::vector< double > shr3d_dEdx = shower3d->dEdx();
      //const int shr3d_bestplane = shower3d->best_plane();

      //           std::cout<<"SHOWER3D_ENERGY: best plane: "<<shr3d_bestplane<<std::endl;
      //for(auto &en:shr3d_energy){
      //    std::cout<<en<<" ";
      //}
      if(shr3d_energy.size()==3){
        vars.m_reco_shower3d_energy_plane0[i_shr] = shr3d_energy[0];
        vars.m_reco_shower3d_energy_plane1[i_shr] = shr3d_energy[1];
        vars.m_reco_shower3d_energy_plane2[i_shr] = shr3d_energy[2];
      }else{
        vars.m_reco_shower3d_energy_plane0[i_shr] =-99;
        vars.m_reco_shower3d_energy_plane1[i_shr] =-99;
        vars.m_reco_shower3d_energy_plane2[i_shr] =-999;
      }

      //         std::cout<<std::endl<<"SHOWER3D_DEDX: "<<std::endl;
      //for(auto &dedx: shr3d_dEdx){
      //    std::cout<<dedx<<" ";
      //}
      if(shr3d_dEdx.size()==3){
        vars.m_reco_shower3d_dEdx_plane0[i_shr] = shr3d_dEdx[0];
        vars.m_reco_shower3d_dEdx_plane1[i_shr] = shr3d_dEdx[1];
        vars.m_reco_shower3d_dEdx_plane2[i_shr] = shr3d_dEdx[2];
      }else{
        vars.m_reco_shower3d_dEdx_plane0[i_shr] =-99;
        vars.m_reco_shower3d_dEdx_plane1[i_shr] =-99;
        vars.m_reco_shower3d_dEdx_plane2[i_shr] =-999;
      }


      //------------- calorimetry ------------
      vars.m_reco_shower_energy_plane0[i_shr] = CalcEShowerPlane(hits, 0, paras);
      vars.m_reco_shower_energy_plane1[i_shr] = CalcEShowerPlane(hits, 1, paras);
      vars.m_reco_shower_energy_plane2[i_shr] = CalcEShowerPlane(hits, 2, paras);

      vars.m_reco_shower_energy_max[i_shr] = std::max( vars.m_reco_shower_energy_plane0[i_shr], std::max( vars.m_reco_shower_energy_plane1[i_shr] , vars.m_reco_shower_energy_plane2[i_shr]));

      vars.m_reco_shower_plane0_nhits[i_shr] = getNHitsPlane(hits, 0);
      vars.m_reco_shower_plane1_nhits[i_shr] = getNHitsPlane(hits, 1);
      vars.m_reco_shower_plane2_nhits[i_shr] = getNHitsPlane(hits, 2);

      vars.m_reco_shower_plane0_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 0);
      vars.m_reco_shower_plane1_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 1);
      vars.m_reco_shower_plane2_meanRMS[i_shr] = getMeanHitWidthPlane(hits, 2);


      //currently only run on 1 shower events
      if(showers.size()==1){
        for(auto &h: hits){ 

          int plane= h->View();
          int wire = h->WireID().Wire;
          int tick = h->PeakTime();

          vars.m_reco_shower_hit_tick.push_back(tick);
          vars.m_reco_shower_hit_plane.push_back(plane);
          vars.m_reco_shower_hit_wire.push_back(wire);
        }
      }


      //std::cout<<"The energy on each plane is 0: "<< vars.m_reco_shower_energy_plane0[i_shr]<<", 1: "<< vars.m_reco_shower_energy_plane1[i_shr]<<", 2: "<<  vars.m_reco_shower_energy_plane2[i_shr]<<std::endl;


      vars.m_reco_shower_dQdx_plane0[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 0 , triggeroffset, theDetector, paras);
      vars.m_reco_shower_dQdx_plane1[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 1 , triggeroffset, theDetector, paras);
      vars.m_reco_shower_dQdx_plane2[i_shr] = CalcdQdxShower(shower,clusters, clusterToHitMap, 2 , triggeroffset, theDetector, paras);
      vars.m_reco_shower_dEdx_plane0[i_shr] = CalcdEdxFromdQdx(vars.m_reco_shower_dQdx_plane0[i_shr], paras);
      vars.m_reco_shower_dEdx_plane1[i_shr] = CalcdEdxFromdQdx(vars.m_reco_shower_dQdx_plane1[i_shr], paras);

      vars.m_reco_shower_dEdx_plane2[i_shr] = CalcdEdxFromdQdx(vars.m_reco_shower_dQdx_plane2[i_shr], paras);

      vars.m_reco_shower_dEdx_plane0_median[i_shr] = getMedian(vars.m_reco_shower_dEdx_plane0[i_shr]);
      vars.m_reco_shower_dEdx_plane1_median[i_shr] = getMedian(vars.m_reco_shower_dEdx_plane1[i_shr]);
      vars.m_reco_shower_dEdx_plane2_median[i_shr] = getMedian(vars.m_reco_shower_dEdx_plane2[i_shr]);


      std::vector< double > reco_shr_angles_wrt_wires;
      for (geo::PlaneGeo const& plane: paras.s_wireReadout->Iterate<geo::PlaneGeo>()) {
        //6 planes in SBND
        //WireAngleToVertical  : 30 ,150,90,150,30 ,90
        //ub wire angles    : 30 ,150,90  (respected to beam,z)
        //Pitch        : 0.3,0.3,0.3,0.3,0.3,0.3

        const double angToVert(paras.s_wireReadout->WireAngleToVertical(plane.View(), plane.ID())+0.5*M_PI);//wire angle respected to z + pi/2

        TVector3 wire_vector;  
        if(abs(angToVert) < 1e-9 ){
          wire_vector = {0,0,1};
        } else{
          wire_vector = { 0 , sin(angToVert) ,cos(angToVert) };
        }

        //    std::cout<<" Angle "<<angToVert<<" Get Vec y="<<wire_vector[1]<< " z= "<<wire_vector[2]<<std::endl;
        reco_shr_angles_wrt_wires.push_back( abs(0.5*M_PI-acos(wire_vector.Dot(shr_dir))) );

        if(reco_shr_angles_wrt_wires.size()==3) break;
      }
      vars.m_reco_shower_angle_wrt_wires_plane0[i_shr] = reco_shr_angles_wrt_wires[0];
      vars.m_reco_shower_angle_wrt_wires_plane1[i_shr] = reco_shr_angles_wrt_wires[1];
      vars.m_reco_shower_angle_wrt_wires_plane2[i_shr] = reco_shr_angles_wrt_wires[2];

      vars.m_reco_shower_dQdx_plane0_median[i_shr] = getMedian(vars.m_reco_shower_dQdx_plane0[i_shr]);
      vars.m_reco_shower_dQdx_plane1_median[i_shr] = getMedian(vars.m_reco_shower_dQdx_plane1[i_shr]);
      vars.m_reco_shower_dQdx_plane2_median[i_shr] = getMedian(vars.m_reco_shower_dQdx_plane2[i_shr]);



      vars.m_reco_shower_dEdx_plane0_mean[i_shr] = std::accumulate(vars.m_reco_shower_dEdx_plane0[i_shr].begin(), vars.m_reco_shower_dEdx_plane0[i_shr].end(), 0.0)/((double)vars.m_reco_shower_dEdx_plane0[i_shr].size()); 
      vars.m_reco_shower_dEdx_plane1_mean[i_shr] = std::accumulate(vars.m_reco_shower_dEdx_plane1[i_shr].begin(), vars.m_reco_shower_dEdx_plane1[i_shr].end(), 0.0)/((double)vars.m_reco_shower_dEdx_plane1[i_shr].size()); 
      vars.m_reco_shower_dEdx_plane2_mean[i_shr] = std::accumulate(vars.m_reco_shower_dEdx_plane2[i_shr].begin(), vars.m_reco_shower_dEdx_plane2[i_shr].end(), 0.0)/((double)vars.m_reco_shower_dEdx_plane2[i_shr].size()); 

      auto maxp0 = std::max_element(vars.m_reco_shower_dEdx_plane0[i_shr].begin(), vars.m_reco_shower_dEdx_plane0[i_shr].end());
      auto maxp1 = std::max_element(vars.m_reco_shower_dEdx_plane1[i_shr].begin(), vars.m_reco_shower_dEdx_plane1[i_shr].end());
      auto maxp2 = std::max_element(vars.m_reco_shower_dEdx_plane2[i_shr].begin(), vars.m_reco_shower_dEdx_plane2[i_shr].end());
      auto minp0 = std::min_element(vars.m_reco_shower_dEdx_plane0[i_shr].begin(), vars.m_reco_shower_dEdx_plane0[i_shr].end());
      auto minp1 = std::min_element(vars.m_reco_shower_dEdx_plane1[i_shr].begin(), vars.m_reco_shower_dEdx_plane1[i_shr].end());
      auto minp2 = std::min_element(vars.m_reco_shower_dEdx_plane2[i_shr].begin(), vars.m_reco_shower_dEdx_plane2[i_shr].end());


      if(maxp0 == vars.m_reco_shower_dEdx_plane0[i_shr].end()){
        vars.m_reco_shower_dEdx_plane0_max[i_shr] = -999; 
      }else{
        vars.m_reco_shower_dEdx_plane0_max[i_shr] = *maxp0; 
      }

      if(maxp1 == vars.m_reco_shower_dEdx_plane1[i_shr].end()){
        vars.m_reco_shower_dEdx_plane1_max[i_shr] = -999; 
      }else{
        vars.m_reco_shower_dEdx_plane1_max[i_shr] = *maxp1; 
      }

      if(maxp2 == vars.m_reco_shower_dEdx_plane2[i_shr].end()){
        vars.m_reco_shower_dEdx_plane2_max[i_shr] = -999; 
      }else{
        vars.m_reco_shower_dEdx_plane2_max[i_shr] = *maxp2; 
      }


      if(minp0 == vars.m_reco_shower_dEdx_plane0[i_shr].end()){
        vars.m_reco_shower_dEdx_plane0_min[i_shr] = -999; 
      }else{
        vars.m_reco_shower_dEdx_plane0_min[i_shr] = *minp0; 
      }

      if(minp1 == vars.m_reco_shower_dEdx_plane1[i_shr].end()){
        vars.m_reco_shower_dEdx_plane1_min[i_shr] = -999; 
      }else{
        vars.m_reco_shower_dEdx_plane1_min[i_shr] = *minp1; 
      }

      if(minp2 == vars.m_reco_shower_dEdx_plane2[i_shr].end()){
        vars.m_reco_shower_dEdx_plane2_min[i_shr] = -999; 
      }else{
        vars.m_reco_shower_dEdx_plane2_min[i_shr] = *minp2; 
      }


      vars.m_reco_shower_dEdx_plane0_nhits[i_shr] = vars.m_reco_shower_dEdx_plane0[i_shr].size();
      vars.m_reco_shower_dEdx_plane1_nhits[i_shr] = vars.m_reco_shower_dEdx_plane1[i_shr].size();
      vars.m_reco_shower_dEdx_plane2_nhits[i_shr] = vars.m_reco_shower_dEdx_plane2[i_shr].size();

      vars.m_reco_shower_dEdx_amalgamated[i_shr] = getAmalgamateddEdx( 
          vars.m_reco_shower_angle_wrt_wires_plane0[i_shr],  
          vars.m_reco_shower_angle_wrt_wires_plane1[i_shr],  
          vars.m_reco_shower_angle_wrt_wires_plane2[i_shr], 
          vars.m_reco_shower_dEdx_plane0_median[i_shr], 
          vars.m_reco_shower_dEdx_plane1_median[i_shr], 
          vars.m_reco_shower_dEdx_plane2_median[i_shr],
          vars.m_reco_shower_dEdx_plane0_nhits[i_shr], 
          vars.m_reco_shower_dEdx_plane1_nhits[i_shr], 
          vars.m_reco_shower_dEdx_plane2_nhits[i_shr] );

      vars.m_reco_shower_dEdx_amalgamated_nhits[i_shr] = getAmalgamateddEdxNHits(
          vars.m_reco_shower_dEdx_amalgamated[i_shr], 
          vars.m_reco_shower_dEdx_plane0_median[i_shr], 
          vars.m_reco_shower_dEdx_plane1_median[i_shr], 
          vars.m_reco_shower_dEdx_plane2_median[i_shr],
          vars.m_reco_shower_dEdx_plane0_nhits[i_shr], 
          vars.m_reco_shower_dEdx_plane1_nhits[i_shr], 
          vars.m_reco_shower_dEdx_plane2_nhits[i_shr] );

      //-------------- Flashes : Was there a flash in the beavars.m_time and if so was it near in Z? --------------------
      double zmin = vars.m_reco_shower_startz[i_shr];
      double zmax = zmin + vars.m_reco_shower_dirz[i_shr]*vars.m_reco_shower_length[i_shr];
      if(zmin > zmax) std::swap(zmin, zmax);

      double ymin = vars.m_reco_shower_starty[i_shr];
      double ymax = zmin + vars.m_reco_shower_diry[i_shr]*vars.m_reco_shower_length[i_shr];
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

      //            if(g_is_verbose) std::cout<<"AnalyzeShowers()\t||\tnumber of flashes: "<< vars.m_reco_num_flashes<< ""<<std::endl;;
      for(int i_flash = 0; i_flash < vars.m_reco_num_flashes; ++i_flash) {

        double const zcenter=vars.m_reco_flash_zcenter[i_flash];
        //                if(g_is_verbose) std::cout<< "AnalyzeShowers()\t||\tflash z center:" <<vars.m_reco_flash_zcenter[i_flash]<< ""<<std::endl;;
        double const ycenter=vars.m_reco_flash_ycenter[i_flash];
        //                if(g_is_verbose) std::cout<< "AnaluzeShowers()\t||\tflash y center:" <<vars.m_reco_flash_ycenter[i_flash]<< ""<<std::endl;;

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
      if(vars.m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_z = -2;
      vars.m_reco_shower_flash_shortest_distz[i_shr]=shortest_dist_to_flash_z;
      vars.m_reco_shower_flash_shortest_index_z[i_shr]=shortest_dist_to_flash_index_z;

      if(vars.m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_y = -2;
      vars.m_reco_shower_flash_shortest_disty[i_shr]=shortest_dist_to_flash_y;
      vars.m_reco_shower_flash_shortest_index_y[i_shr]=shortest_dist_to_flash_index_y;
      vars.m_reco_shower_flash_shortest_distyz[i_shr]=shortest_dist_to_flash_yz;
      vars.m_reco_shower_flash_shortest_index_yz[i_shr]=shortest_dist_to_flash_index_yz;
      if(vars.m_reco_num_flashes_in_beamgate == 0) shortest_dist_to_flash_yz = -2;

      //end optical flash code


      vars.m_reco_shower_num_daughters[i_shr] = pfp->NumDaughters();  //corresponding PFParticle
      //      std::cout<<" CHECK numebr "<<vars.m_reco_shower_num_daughters[i_shr]<<std::endl;
      if(vars.m_reco_shower_num_daughters[i_shr]>0){
        //currently just look at 1 daughter
        //vars.m_reco_shower_daughter_trackscore[i_shr] = PFPToTrackScoreMap[pfParticleMap[pfp->Daughters().front()]];
        int pfp_size = all_PPFPs.size();
        for(int index = 0; index < pfp_size; index++){
          //          std::cout<<"CHECK Compare "<<pfp->Daughters().front()<<
          //          " "<<all_PPFPs[index].pPFParticle->Self()<<std::endl;
          if( (pfp->Daughters().front()) == all_PPFPs[index].pPFParticle->Self()){
            vars.m_reco_shower_daughter_trackscore[i_shr] = all_PPFPs[index].get_TrackScore();
            break;
          }
        }
      }


      //------------and finally some slice info-----------------

      vars.m_reco_shower_sliceId[i_shr] = ppfp->get_SliceID();//PFPToSliceIdMap[pfp];
      vars.m_reco_shower_nuscore[i_shr] = ppfp->get_NuScore();//sliceIdToNuScoreMap[ vars.m_reco_shower_sliceId[i_shr]] ;
      vars.m_reco_shower_isclearcosmic[i_shr] = ppfp->get_IsClearCosmic();//PFPToClearCosmicMap[pfp];
      vars.m_reco_shower_is_nuslice[i_shr] = ppfp->get_IsNuSlice();//PFPToNuSliceMap[pfp];

      vars.m_reco_shower_trackscore[i_shr] = ppfp->get_TrackScore();
      vars.m_reco_shower_pfparticle_pdg[i_shr] = ppfp->get_PdgCode();

      //            if ( vars.m_reco_shower_sliceId[i_shr] >0) std::cout<<"AnalyzeShowers()\t||\t On Shower: "<<i_shr<<". Pfp id = "<< pfp->Self()<<". The slice id for this shower is "<< vars.m_reco_shower_sliceId[i_shr]<<", the neutrino score for this slice is "<< vars.m_reco_shower_nuscore[i_shr]<<", and is_nuslice = "<<  vars.m_reco_shower_is_nuslice[i_shr]<<". The track score is : "<< vars.m_reco_shower_trackscore[i_shr]<<std::endl;

      i_shr++;

      //std::vector<int> spacers = Printer_header({"Slice","pfpID","Start(x,  ","   y,      ",",      z  )", "trackScore", "Pdg"});
      Printer_content(
          {std::to_string(ppfp->get_SliceID()),
          std::to_string(ppfp->get_PFParticleID()),
          std::to_string(shr_start.X()),
          std::to_string(shr_start.Y()),
          std::to_string(shr_start.Z()),
          std::to_string(ppfp->get_TrackScore()),
          std::to_string(ppfp->get_PdgCode())
          },spacers);

    }

    //Lets sort and order the showers
    vars.m_reco_shower_ordered_energy_index = sort_indexes(vars.m_reco_shower_energy_max);
  }
}
