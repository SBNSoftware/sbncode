#include "sbncode/SinglePhotonAnalysis/SEAview/SEAviewer.h"

#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include <algorithm>
#include <cmath>

namespace seaview{

  int cluster::InNuSlice(const std::map<int, std::vector<art::Ptr<recob::Hit>> >& sliceIDToHitsMap, int nuSliceID){

    const std::vector<art::Ptr<recob::Hit>>& slice_hits = sliceIDToHitsMap.at(nuSliceID);
    bool found_hit_in_slice = false, found_hit_not_in_slice = false;
    for(auto hit : f_hits){
      auto iter = std::find(slice_hits.begin(), slice_hits.end(), hit);
      if( iter  == slice_hits.end())
        found_hit_not_in_slice = true;
      else{
        found_hit_in_slice = true;
      }

      if(found_hit_in_slice && found_hit_not_in_slice) return 0;
    }

    return found_hit_in_slice ? 1 : -1;
  }


  void SEAviewer::TrackLikeClusterAnalyzer(cluster &cl, const std::vector<double> &shower_start_pt_2D, const std::vector<double> &shower_other_pt_2D){

    //first round, grab min ioc, impact, conversion distance, and ADC_histogram, and indx of the hits with min/max IOC
    double max_ioc_to_shower_start = DBL_MIN;
    size_t num_hits = cl.f_hits.size();


    for(size_t i = 0; i!= num_hits; ++i){
      auto &h = cl.f_hits[i]; //type of h: art_ptr<recob::Hit>
      double h_tick = (double)h->PeakTime();
      double h_wire = (double)h->WireID().Wire;

      //geometric properties
      double impact_parameter_to_shower = dist_line_point(shower_start_pt_2D, shower_other_pt_2D, {h_wire, h_tick});
      double conversion_dist_to_shower_start = dist_point_point(h_wire, h_tick, shower_start_pt_2D[0], shower_start_pt_2D[1]);
      double ioc_to_shower_start = impact_parameter_to_shower/conversion_dist_to_shower_start;
      cl.f_min_impact_parameter_to_shower = std::min(cl.f_min_impact_parameter_to_shower, impact_parameter_to_shower);
      cl.f_min_conversion_dist_to_shower_start = std::min(cl.f_min_conversion_dist_to_shower_start, conversion_dist_to_shower_start);

      //remember two hits with min/max IOC
      if( ioc_to_shower_start < cl.f_min_ioc_to_shower_start){
        cl.f_min_ioc_to_shower_start = ioc_to_shower_start;
        cl.start_hit_idx = i; 
      }
      if( ioc_to_shower_start > max_ioc_to_shower_start){ //be careful, should not be "else if" here.
        max_ioc_to_shower_start = ioc_to_shower_start;
        cl.end_hit_idx = i;
      }


    } //end of hit loop


    // second round: group hits in two categories, first half and second half
    // get the direction of the cluster
    auto start_hit_ptr = cl.f_hits.at(cl.start_hit_idx), end_hit_ptr = cl.f_hits.at(cl.end_hit_idx);
    std::vector<double> start_hit_point = { (double)start_hit_ptr->WireID().Wire, (double)start_hit_ptr->PeakTime()};  
    std::vector<double> end_hit_point= { (double)end_hit_ptr->WireID().Wire, (double)end_hit_ptr->PeakTime()};  
    std::vector<double> mid_point = { (start_hit_point[0] + end_hit_point[0])/2, (start_hit_point[1] + end_hit_point[1])/2};
    std::vector<double> start_to_mid_vec = { (- start_hit_point[0] + end_hit_point[0])/2, ( - start_hit_point[1] + end_hit_point[1])/2};
    cl.f_ioc_based_length = sqrt(pow((start_hit_point[0]-end_hit_point[0])*wire_con, 2.0) + pow((start_hit_point[1]-end_hit_point[1])*tick_con, 2.0));

    cl.f_hit_group.resize(num_hits);  
    for(size_t i = 0; i!=num_hits; ++i){
      auto h = cl.f_hits[i]; //type of h: art_ptr<recob::Hit>
      double h_tick = (double)h->PeakTime();
      double h_wire = (double)h->WireID().Wire;

      std::vector<double> mid_to_h_vec = { h_wire - mid_point[0], h_tick - mid_point[1]};

      if( start_to_mid_vec[0]*mid_to_h_vec[0]*pow(wire_con, 2.0) + start_to_mid_vec[1] * mid_to_h_vec[1] *pow(tick_con, 2.0)<= 0 ){
        cl.f_hit_group[i] = 1;
        cl.f_mean_ADC_first_half += h->ROISummedADC();
      }
      else{
        cl.f_hit_group[i] = 2;
        cl.f_mean_ADC_second_half += h->ROISummedADC();
      }
    }
    cl.f_track_treated = true;

    size_t nhits_first_half = std::count(cl.f_hit_group.begin(), cl.f_hit_group.end(), 1);
    size_t nhits_second_half = std::count(cl.f_hit_group.begin(), cl.f_hit_group.end(), 2);
    if(nhits_first_half) cl.f_mean_ADC_first_half /= nhits_first_half;
    if(nhits_second_half) cl.f_mean_ADC_second_half /= nhits_second_half;
    cl.f_mean_ADC_first_to_second_ratio = (cl.f_mean_ADC_second_half != 0 ? cl.f_mean_ADC_first_half/cl.f_mean_ADC_second_half : 0.0);

    if( num_hits < 2){
      std::cerr << "SEAviewer::TrackLikeClusterAnalyzer\t|| this cluster only has " << num_hits << " hits, can't calculate direction, skipping it... " << std::endl;
      return;
    }

    //angle between the track cluster and the shower direction
    //cluster direction unit vector
    double start_to_mid_length = sqrt( pow(start_to_mid_vec[0]*wire_con, 2.0) + pow(start_to_mid_vec[1]*tick_con, 2.0));
    std::vector<double> cluster_dir = {start_to_mid_vec[0]*wire_con/start_to_mid_length, start_to_mid_vec[1]*tick_con/start_to_mid_length};
    //shower direction unit vector
    double shower_direction_length = sqrt(pow((shower_start_pt_2D[0] - shower_other_pt_2D[0])*wire_con, 2.0) + pow((shower_start_pt_2D[1] - shower_other_pt_2D[1])*tick_con, 2.0));
    std::vector<double> shower_dir = { (shower_other_pt_2D[0] - shower_start_pt_2D[0])*wire_con/shower_direction_length, (shower_other_pt_2D[1] - shower_start_pt_2D[1])*tick_con/shower_direction_length};
    //angle between two unit vector, in radian
    cl.f_angle_wrt_shower_direction = acos( cluster_dir[0]*shower_dir[0] + cluster_dir[1]*shower_dir[1]);


    if(cl.f_score.min_wire == cl.f_score.max_wire){
      std::cout << "SEAviewer::TrackLikeClusterAnalyzer\t|| this cluster spans only on 1 wire: " << cl.f_score.min_wire << ", setting its straight-line fit chi2 to 0.." << std::endl;
      return;
    }
    //fit to wire-tick plot of the cluster, see how straight cluster is
    TF1 *f1 = new TF1("f1", "1 ++ x", cl.f_score.min_wire, cl.f_score.max_wire);
    //TGraph* graph_copy = (TGraph*)cl.f_graph.Clone("temp");
    //int fit_status = graph_copy->Fit(f1, "RQ0"); //if fit status is 0, the fit is ok.
    int fit_status = cl.f_graph.Fit(f1, "RQ0"); //if fit status is 0, the fit is ok.
    if(fit_status == 0){
      //f1 = graph_copy->GetFunction("f1");
      f1 = cl.f_graph.GetFunction("f1");
      cl.f_fit_chi2 = f1->GetChisquare();
    }
    //std::cout << "SEAviewer::TrackLikeClusterAnalyzer\t|| End" << std::endl;
  }


  // constructor
  // CHECK, this constructor is not found;
  SEAviewer::SEAviewer(std::string intag, 
      geo::GeometryCore const * ingeom, 
      geo::WireReadoutGeom const * inwireReadout,
      detinfo::DetectorPropertiesData const & intheDetector )
    : tag(intag), geom(ingeom), wireReadout{inwireReadout}, theDetector(intheDetector){
      chan_max = {-9999,-9999,-9999};
      chan_min = {9999,9999,9999};
      tick_max = -99999;
      tick_min = 99999;

      plot_true_vertex = false;
      vertex_tick.resize(3); 
      vertex_chan.resize(3); 
      vertex_graph.resize(3); 

      true_vertex_tick.resize(3); 
      true_vertex_chan.resize(3); 
      true_vertex_graph.resize(3); 

      tick_shift = 350;
      chan_shift = 100;

      n_showers=0;
      n_pfps = 0;
      has_been_clustered = false;
      hit_threshold = -10;

      rangen = new TRandom3(0); //same seed everytime
    }


  int SEAviewer::setBadChannelList(std::vector<std::pair<int,int>> &in){
    m_bad_channel_list = in;
    return 0;
  }

  int SEAviewer::addHitsToConsider(std::vector<art::Ptr<recob::Hit>>& hits){
    for(auto &h: hits){
      map_unassociated_hits[h] = true;
      map_considered_hits[h] = true;
    }
    return 0;
  }

  int SEAviewer::filterConsideredHits(double dist_to_vertex){

    //collect all hits that are under consideration for clustering
    std::vector<art::Ptr<recob::Hit>> current_hits;
    for(auto map_iter = map_considered_hits.begin(); map_iter != map_considered_hits.end(); ++map_iter){
      current_hits.push_back(map_iter->first);
    }

    //remove hits that are too far from vertex
    for(auto &h: current_hits){
      int p = h->View();
      double wire = (double)h->WireID().Wire;
      double tick = (double)h->PeakTime();
      double dist = dist_point_point(wire, tick, vertex_chan[p], vertex_tick[p]);

      if(dist > dist_to_vertex){
        map_unassociated_hits.erase(h);
        map_considered_hits.erase(h);
      }
    }
    return 0;
  }

  int SEAviewer::setHitThreshold(double h){
    hit_threshold = h;
    return 0;    
  }

  int SEAviewer::addAllHits(std::vector<art::Ptr<recob::Hit>>& hits){

    std::vector<std::vector<double>>  vec_tick(3);
    std::vector<std::vector<double>>  vec_chan(3);

    for(auto&h: hits){
      if(map_considered_hits.count(h)==0){   // if h is not in the map, push its plane ID, wire ID, and time tick to the vectors
        double wire = (double)h->WireID().Wire;
        vec_chan[(int)h->View()].push_back(wire);
        vec_tick[(int)h->View()].push_back((double)h->PeakTime());
        //tick_max = std::max(tick_max, (double)h->PeakTime());
        //tick_min = std::min(tick_min, (double)h->PeakTime());
        //chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
        //chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);
      }
    }

    for(int i=0; i<3; i++){
      vec_all_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); //implicitly converted to TGraph
    }

    vec_all_ticks = vec_tick;
    vec_all_chans = vec_chan;

    return 0;
  }

  std::vector<int> SEAviewer::calcUnassociatedHits(){
    int n_unassoc=0;
    int n_below_thresh = 0;
    std::vector<std::vector<double>>  vec_tick(3);
    std::vector<std::vector<double>>  vec_chan(3);
    std::vector<std::vector<std::vector<double>>> vec_pts(3);
    std::vector<std::vector<art::Ptr<recob::Hit>>> vec_hits(3);


    int n_all =map_considered_hits.size();

    for(const auto &pair: map_considered_hits ){
      auto& h = pair.first;  //type of h: recob::Hit
      if(map_unassociated_hits.count(h) !=0 && map_unassociated_hits[h]){

        if(h->ROISummedADC() < hit_threshold){
          n_below_thresh++;
          continue;
        }

        // if summed ADC of the hit passes threshold                
        n_unassoc++;
        double wire = (double)h->WireID().Wire;
        double tick = (double)h->PeakTime();
        vec_chan[(int)h->View()].push_back(wire);
        vec_tick[(int)h->View()].push_back(tick);

        vec_pts[(int)h->View()].push_back({wire,tick});
        vec_hits[(int)h->View()].push_back(h);
        tick_max = std::max(tick_max, tick);
        tick_min = std::min(tick_min, tick);
        //chan_max, chan_min stores: max, min unassociated channel for each plane
        chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
        chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

      }

    }

    for(int i=0; i<3; i++){
      vec_unass_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); 
    }

    vec_unass_ticks = vec_tick;
    vec_unass_chans = vec_chan;
    vec_unass_pts = vec_pts;
    vec_unass_hits = vec_hits;

    return {n_all,n_unassoc,n_below_thresh};
  }


  int SEAviewer::addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits, std::string legend, double arg1, double arg2, double arg3){
    n_pfps++;

    format_legend(legend, arg1, arg2, arg3);

    vec_pfp_legend.push_back(legend);

    std::vector<std::vector<double>>  vec_tick(3);
    std::vector<std::vector<double>>  vec_chan(3);

    for(auto &h: hits){
      double wire = (double)h->WireID().Wire;
      vec_chan[(int)h->View()].push_back(wire);
      vec_tick[(int)h->View()].push_back((double)h->PeakTime());
      tick_max = std::max(tick_max, (double)h->PeakTime());
      tick_min = std::min(tick_min, (double)h->PeakTime());
      chan_max[(int)h->View()] = std::max( chan_max[(int)h->View()],wire);
      chan_min[(int)h->View()] = std::min( chan_min[(int)h->View()],wire);

      //remove from unassociated hits
      map_unassociated_hits[h] = false;
    }

    std::vector<TGraph> t_graphs;
    for(int i=0; i<3; i++){
      t_graphs.emplace_back(vec_tick[i].size(),&vec_chan[i][0],&vec_tick[i][0]); 
    }
    vec_graphs.push_back(t_graphs);
    vec_ticks.push_back(vec_tick);
    vec_chans.push_back(vec_chan);
    return 0;
  }

  int SEAviewer::addShower(art::Ptr<recob::Shower>&shr){
    n_showers++;

    vec_showers.push_back(shr); 

    return 0;
  }

  int SEAviewer::addTrack(art::Ptr<recob::Track>&trk){
    n_tracks++;

    vec_tracks.push_back(trk); 

    return 0;
  }


  std::vector<std::vector<double>> SEAviewer::to2D(std::vector<double> & threeD){

    auto const& TPC = *geom->begin<geo::TPCGeo>();  // first TPC of detector
    auto ID = TPC.ID(); 
    int fCryostat = ID.Cryostat;
    int fTPC = ID.TPC;

    std::vector<std::vector<double>> ans(3);

    for(int i=0; i<3; i++){
      double wire = (double)calcWire(threeD[1], threeD[2], i, fTPC, fCryostat, *wireReadout);
      double time = calcTime(threeD[0], i, fTPC,fCryostat, theDetector);

      ans[i] = {wire,time};
    }

    return ans;
  }



  int SEAviewer::loadVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z){

    auto const& TPC = *geom->begin<geo::TPCGeo>();
    auto ID = TPC.ID();
    int fCryostat = ID.Cryostat;
    int fTPC = ID.TPC;

    for(int i=0; i<3; i++){

      // use vector here, so that to plot the single point using TGraph
      std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *wireReadout)};
      std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, theDetector)};

      vertex_tick[i] = time[0];
      vertex_chan[i] = wire[0];

      chan_max[i] = std::max( chan_max[i],wire[0]);
      chan_min[i] = std::min( chan_min[i],wire[0]);

      TGraph gtmp(1, &wire[0], &time[0]); 
      vertex_graph[i] = gtmp;
    }

    return 0;
  }


  int SEAviewer::addTrueVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z){

    plot_true_vertex = true;

    auto const& TPC = *geom->begin<geo::TPCGeo>();
    auto ID = TPC.ID();
    int fCryostat = ID.Cryostat;
    int fTPC = ID.TPC;

    for(int i=0; i<3; i++){

      std::vector<double> wire = {(double)calcWire(m_vertex_pos_y, m_vertex_pos_z, i, fTPC, fCryostat, *wireReadout)};
      std::vector<double> time = {calcTime(m_vertex_pos_x, i, fTPC,fCryostat, theDetector)};

      true_vertex_tick[i] = time[0];
      true_vertex_chan[i] = wire[0];

      chan_max[i] = std::max( chan_max[i],wire[0]);
      chan_min[i] = std::min( chan_min[i],wire[0]);

      TGraph gtmp(1, &wire[0], &time[0]); 
      true_vertex_graph[i] = gtmp;
    }

    return 0;
  }


  int SEAviewer::Print(double plot_distance){


    std::string print_name = "SEAview_"+tag;
    TCanvas *can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,800);
    can->Divide(4,1,0,0.1);

    double plot_point_size=0.4;        

    //******************************* First plot "Vertex" ***************************************

    //Calculate some things
    //Guanqun: what does tick_min - tick_shift actually mean?
    double real_tick_min =  (fabs(vertex_tick[0] - (tick_min-tick_shift))*tick_con > plot_distance)  ? vertex_tick[0]-plot_distance/tick_con  : tick_min-tick_shift  ;
    double real_tick_max =  (fabs(vertex_tick[0] - (tick_max+tick_shift))*tick_con > plot_distance)  ? vertex_tick[0]+plot_distance/tick_con  : tick_max+tick_shift  ;
    //double real_tick_min = tick_min-tick_shift  ;
    //double real_tick_max = tick_max+tick_shift  ;


    std::vector<double> real_wire_min(3); //real x axis edges for 3 planes
    std::vector<double> real_wire_max(3);

    for(int i=0; i<3; i++){
      TPad * pader = (TPad*)can->cd(i+1);

      if(i==0 || i ==4 || i == 8) pader->SetLeftMargin(0.1);

      //only show area surrounding the vertex up to std::min(plot_distance, distance_bw_vertex_channel_min/max)
      real_wire_min[i] =  (fabs(vertex_chan[i] - (chan_min[i]-chan_shift))*wire_con > plot_distance ) ? vertex_chan[i]-plot_distance/wire_con  : chan_min[i]-chan_shift  ;
      real_wire_max[i] =  (fabs(vertex_chan[i] - (chan_max[i]+chan_shift))*wire_con > plot_distance ) ? vertex_chan[i]+plot_distance/wire_con  : chan_max[i]+chan_shift  ;

      //fix the area to show, always show area large enough to hold all track/showers
      //real_wire_min[i] =   chan_min[i]-chan_shift  ;
      //real_wire_max[i] =   chan_max[i]+chan_shift  ;

      vertex_graph[i].SetMarkerStyle(29);
      vertex_graph[i].SetMarkerSize(2);
      vertex_graph[i].SetMarkerColor(kMagenta-3);
      vertex_graph[i].GetYaxis()->SetRangeUser(real_tick_min,real_tick_max);
      vertex_graph[i].GetXaxis()->SetLimits(real_wire_min[i], real_wire_max[i]);
      vertex_graph[i].SetTitle(("Plane " +std::to_string(i)).c_str());
      vertex_graph[i].GetYaxis()->SetTitle("Peak Hit Time Tick");
      vertex_graph[i].GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
      vertex_graph[i].Draw("ap");

      if(i>0){
        vertex_graph[i].GetYaxis()->SetLabelOffset(999);
        vertex_graph[i].GetYaxis()->SetLabelSize(0);
      }
    }

    /********************************* Non Slice  Hits ****************************/


    for(int i=0; i<3; i++){
      can->cd(i+1);
      if(vec_all_graphs[i].GetN()>0){//need a check in case this track has no hits on this plane.
        vec_all_graphs[i].Draw("p same"); 
        vec_all_graphs[i].SetMarkerColor(kGray);
        vec_all_graphs[i].SetFillColor(kWhite);
        vec_all_graphs[i].SetMarkerStyle(20);
        vec_all_graphs[i].SetMarkerSize(plot_point_size*0.75);
      }
    }

    //******************************** DeadWireRegions********************************************
    for(size_t i=0; i< m_bad_channel_list.size(); i++){
      int badchan = m_bad_channel_list[i].first;                                       
      int ok = m_bad_channel_list[i].second;       

      if(ok>1)continue;
      auto hs = wireReadout->ChannelToWire(badchan); //type of hs: vector containing the ID of all the connected wires

      int thisp = (int)hs[0].Plane;
      double bc = hs[0].Wire;

      if(real_wire_min[thisp] < bc && bc < real_wire_max[thisp] ){
        //if(chan_min[thisp]-chan_shift < bc && bc < chan_max[thisp]+chan_shift ){}
        can->cd(thisp+1);
        TLine *l = new TLine(bc, real_tick_min, bc, real_tick_max);
        //TLine *l = new TLine(bc,tick_min-tick_shift,bc,tick_max+tick_shift);
        l->SetLineColor(kGray+1);
        l->Draw("same");
        //can->cd(thisp+5);// Guanqun: how many values can plane ID take?
        //l->Draw("same");
        //can->cd(thisp+9);
        //l->Draw("same");
      }
    }


    ///******************************** Plotting all PFP's *********************************8

    std::vector<int> tcols = {kRed-7, kBlue-7, kGreen-3, kOrange-3, kCyan-3, kMagenta-3, kGreen+1 , kRed+1};
    int used_col=0;

    if(n_pfps > (int)tcols.size()){
      for(int i =0; i< (int)(n_pfps +2); i++){
        //tcols.push_back(tcols[(int)rangen->Uniform(0,7)]+(int)rangen->Uniform(-5,5));
        tcols.push_back(kRed);
      }
    }


    for(int p=0; p<n_pfps; p++){

      int tcol = tcols[used_col];
      used_col++;

      for(int i=0; i<3; i++){
        can->cd(i+1);
        if(vec_graphs[p][i].GetN()>0){//need a check in case this track has no hits on this plane.

          vec_graphs[p][i].Draw("p same"); 
          vec_graphs[p][i].SetMarkerColor(tcol);
          vec_graphs[p][i].SetFillColor(tcol);
          vec_graphs[p][i].SetMarkerStyle(20);
          vec_graphs[p][i].SetMarkerSize(plot_point_size);
        }
      }
    }


    //Plot all Shower lines. Might need a bit of work here..
    std::vector<TLine*> lines;  

    for(size_t s=0; s<vec_showers.size(); ++s){
      std::vector<double> shr_start_3D= {vec_showers[s]->ShowerStart().X(), vec_showers[s]->ShowerStart().Y(),vec_showers[s]->ShowerStart().Z()};
      std::vector<double> shr_other_3D=  {vec_showers[s]->ShowerStart().X()+vec_showers[s]->Direction().X(),vec_showers[s]->ShowerStart().Y()+vec_showers[s]->Direction().Y(), vec_showers[s]->ShowerStart().Z()+vec_showers[s]->Direction().Z()};

      //std::cout<<" "<<shr_start_3D[0]<<" "<<shr_start_3D[1]<<" "<<shr_start_3D[2]<<std::endl;
      //std::cout<<" "<<shr_other_3D[0]<<" "<<shr_other_3D[1]<<" "<<shr_other_3D[2]<<std::endl;

      std::vector<std::vector<double>> start_pt =   this->to2D(shr_start_3D);
      std::vector<std::vector<double>> other_pt =   this->to2D(shr_other_3D);

      for(int i=0; i<3; i++){
        //std::cout<<start_pt[i][0]<<" "<<start_pt[i][1]<<" "<<other_pt[i][0]<<" "<<other_pt[i][1]<<std::endl;
        can->cd(i+1);
        double slope = (start_pt[i][1]-other_pt[i][1])/(start_pt[i][0]-other_pt[i][0]);
        double inter = start_pt[i][1]-slope*start_pt[i][0];

        double x1_plot = other_pt[i][0];//chan_min[i]-chan_shift;
        double y1_plot = slope*x1_plot+inter;

        double x2_plot;
        if(other_pt[i][0]<start_pt[i][0]){
          //x2_plot = chan_max[i]+chan_shift; //guanqun: my guess is this needs to be updated as well to use real_wire_max/min
          x2_plot = real_wire_max[i];
        }else{
          //x2_plot = chan_min[i]-chan_shift;
          x2_plot = real_wire_min[i];    
        }
        double y2_plot = slope*x2_plot+inter;

        can->cd(i+1);
        TLine *l = new TLine(x1_plot, y1_plot, x2_plot, y2_plot);
        lines.push_back(l);
        l->SetLineColorAlpha(tcols[s],0.5);
        l->SetLineWidth(1);
        l->SetLineStyle(2);
        l->Draw();

      }

    }

    /********************************* Unassociated Hits ****************************/
    for(int i=0; i<3; i++){
      can->cd(i+1);
      if(vec_unass_graphs[i].GetN()>0){//need a check in case this track has no hits on this plane.

        vec_unass_graphs[i].Draw("p same"); 
        vec_unass_graphs[i].SetMarkerColor(kBlack);
        vec_unass_graphs[i].SetFillColor(kBlack);
        vec_unass_graphs[i].SetMarkerStyle(20);
        vec_unass_graphs[i].SetMarkerSize(plot_point_size);
      }
    }

    /******************************* Clustered Hits ***********************************/
    // draw cluster hits after drawing all unassociated hits such that clustered hits would be colored while un-clustered ones will be black.
    if(has_been_clustered){ 

      std::vector<int> cluster_colors(vec_clusters.size()+1,0);
      std::vector<int> base_col = {632,416, 600, 400, 616,  432,  800, 820,  840, 860,  880, 900};

      for(size_t j=0; j< vec_clusters.size()+1; j++){
        int b = (int)rangen->Uniform(0,11);
        int mod = (int)rangen->Uniform(-10,+3);

        cluster_colors[j] = base_col[b]+mod;
      }
      int c_offset = 0;

      for(auto &c: vec_clusters){
        int pl = c.getPlane();
        can->cd(pl+1);   
        if (c.getGraph()->GetN()>0){
          c.getGraph()->Draw("p same");
          c.getGraph()->SetMarkerColor(cluster_colors[c_offset]);
          c.getGraph()->SetFillColor(cluster_colors[c_offset]);
          c.getGraph()->SetMarkerStyle(20);
          //c.getGraph()->SetMarkerSize(plot_point_size);
          c.getGraph()->SetMarkerSize(plot_point_size*1.5);
          //std::cout<<"Printing cluster "<<c.getID()<<" on plane "<<pl<<" col "<<cluster_colors[c_offset]<<std::endl;
          //auto ll = c.getPTS();
          //for(auto &p :ll){
          //    std::cout<<p[0]<<" "<<p[1]<<std::endl;
          //}
        }
        c_offset++;
      }
    }//end clusters



    //****** just plto vertex again with elipse;
    for(int i=0; i<3; i++){
      can->cd(i+1);
      vertex_graph[i].Draw("p same");

      double rad_cm = 12.0;
      TEllipse ell_p(vertex_chan[i],vertex_tick[i],rad_cm/wire_con,rad_cm/tick_con);
      ell_p.SetLineColor(kRed);
      ell_p.SetFillStyle(0);
      ell_p.Draw("same");
    }

    //**************************** INFO ***************************/
    TPad *p_top_info = (TPad*)can->cd(4);
    p_top_info->cd();

    /*TLatex pottex;
      pottex.SetTextSize(0.045);
      pottex.SetTextAlign(13);  //align at top
      pottex.SetNDC();
      std::string pot_draw = "Run: "+std::to_string(m_run_number)+" SubRun: "+std::to_string(m_subrun_number)+" Event: "+std::to_string(m_event_number);
      pottex.DrawLatex(.1,.94, pot_draw.c_str());
      */
    TLegend l_top(0.1,0.0,0.9,1.0);
    l_top.SetTextSize(0.05);

    for(int p=0; p<n_pfps; p++){


      if(vec_graphs[p][0].GetN()>0){
        l_top.AddEntry(&vec_graphs[p][0],vec_pfp_legend[p].c_str(),"f");
      }else if(vec_graphs[p][1].GetN()>0){
        l_top.AddEntry(&vec_graphs[p][1],vec_pfp_legend[p].c_str(),"f");
      }else if(vec_graphs[p][2].GetN()>0){
        l_top.AddEntry(&vec_graphs[p][2],vec_pfp_legend[p].c_str(),"f");
      }

    }

    // draw legend for clustered hits if there is any
    for(const auto &cluster : vec_clusters){

      // only consider clusters that are second shower candidates
      if(cluster.getLegend().empty()) continue;

      // if the cluster is out of the plotting range, do not include it in the legend
      if(cluster.InRange(real_tick_max, real_tick_min, real_wire_max[cluster.getPlane()], real_wire_min[cluster.getPlane()])){
        l_top.AddEntry(cluster.getGraph(), cluster.getLegend().c_str(), "f");
      }
    }

    l_top.SetHeader(print_name.c_str(),"C");
    l_top.SetLineWidth(0);
    l_top.SetLineColor(kWhite);
    l_top.Draw("same");



    can->Update();
    can->SaveAs((print_name+".pdf").c_str(),"pdf");


    return 0;
  }

  int SEAviewer::runseaDBSCAN(double min_pts, double eps){

    has_been_clustered = true;
    num_clusters = {0,0,0};
    cluster_labels.resize(3);

    for(int i=0; i<3; i++){

      std::cout<<"SinglePhoton::seaDBSCAN\t||\tStarting to run seaDBSCAN for plane: "<<i<<" has "<<vec_unass_pts[i].size()<<" pts to do using eps: "<<eps<<" and min_pts: "<<min_pts<<std::endl; 
      seaDBSCAN ReCluster(eps,min_pts);
      cluster_labels[i] =  ReCluster.Scan2D(vec_unass_pts[i]);

      for(auto &c: cluster_labels[i]){
        num_clusters[i] = std::max(c,num_clusters[i]);

      }

      std::cout << "SinglePhoton::seaDBSCAN\t||\tOn this plane "<<i<<" seaDBSCAN found: "<<num_clusters[i]<<" clusters"<<std::endl;
    }

    //Now fill the clusters

    for(int i=0; i<3; i++){

      for(int c=0; c<num_clusters[i]+1; c++){

        std::vector<std::vector<double>> t_pts;
        std::vector<art::Ptr<recob::Hit>> hitz;


        for(size_t p=0; p< vec_unass_pts[i].size(); p++){
          if(cluster_labels[i][p] == 0) continue;//noise 
          if(cluster_labels[i][p] == c){

            t_pts.push_back(vec_unass_pts[i][p]);
            hitz.push_back(vec_unass_hits[i][p]);
          }

        }

        if(hitz.size()!=0){
          std::cout<<"SinglePhoton::seaDBSCAN\t||\t Cluster "<<vec_clusters.size()+1<<" has "<<hitz.size()<<" hitz on plane "<<i<<std::endl;
          vec_clusters.emplace_back(vec_clusters.size()+1,i,t_pts,hitz);
        }
      }
    }

    return 0;
  }

  std::vector<double> SEAviewer::analyzeTrackLikeClusters(double eps, const std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap, const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,std::vector<seaview::cluster>& vec_in_clusters ){


    // Grab the shower start and shower direction
    std::vector<double> shr_start_3D= {vec_showers[0]->ShowerStart().X(), vec_showers[0]->ShowerStart().Y(),vec_showers[0]->ShowerStart().Z()};
    std::vector<double> shr_other_3D=  {vec_showers[0]->ShowerStart().X()+vec_showers[0]->Direction().X(),vec_showers[0]->ShowerStart().Y()+vec_showers[0]->Direction().Y(), vec_showers[0]->ShowerStart().Z()+vec_showers[0]->Direction().Z()};

    std::vector<std::vector<double>> shr_start_pt =   this->to2D(shr_start_3D);
    std::vector<std::vector<double>> shr_other_pt =   this->to2D(shr_other_3D);



    //Loop over all clusters
    for(size_t c=0; c< vec_clusters.size(); c++){

      auto hitz = vec_clusters[c].getHits(); // type of hitz: std::vector<art::Ptr<recob::Hit>>
      int num_hits_in_cluster = vec_clusters[c].getHits().size();
      int pl = vec_clusters[c].getPlane();

      //Need to modify this a bit
      auto ssscorz = SeaviewScoreCluster(pl,c+1, hitz ,vertex_chan[pl], vertex_tick[pl], vec_showers.front());
      vec_clusters[c].setScore(ssscorz);
      vec_clusters[c].setWireTickBasedLength(sqrt(pow((ssscorz.max_wire-ssscorz.min_wire)*wire_con, 2.0) + pow((ssscorz.max_tick - ssscorz.min_tick)*tick_con, 2.0)));
      BasicClusterCalorimetry(vec_clusters[c]);
      TrackLikeClusterAnalyzer(vec_clusters[c], shr_start_pt[pl], shr_other_pt[pl]);

      //This is just checking if its in, we can do this earlier; TODO
      //TODO: is_in_shower, get back to primary shower (at least available)
      //Sim Stuff
      //Draw Direction on plot
      //Delauney on here might be good, that said, we have a LOT of things. Hmm, cap at 50 hits maybe? 
      int is_in_shower = SeaviewCompareToShowers(pl,c+1, hitz ,vertex_chan[pl], vertex_tick[pl], vec_showers, showerToPFParticleMap, pfParticleToHitsMap,eps);
      vec_clusters[c].setShowerRemerge(is_in_shower);

      std::string sname = "Cluster "+std::to_string(c)+", Hits: "+std::to_string(num_hits_in_cluster)+", PCA "+std::to_string(ssscorz.pca_0)+", Theta:" +std::to_string(ssscorz.pca_theta)+", Wires: "+std::to_string(ssscorz.n_wires)+ ", Ticks: "+std::to_string(ssscorz.n_ticks)+", ReMerged: "+std::to_string(is_in_shower);
      std::cout<<sname << "\n" <<std::endl;


    }//cluster loop

    vec_in_clusters = vec_clusters;

    return {0};

  }

  std::vector<double> SEAviewer::analyzeShowerLikeClusters(double eps, const std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap, const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,std::vector<seaview::cluster>& vec_in_clusters ){

    /*
       std::vector<std::vector<double>> percent_matched(vec_clusters.size(),  std::vector<double>(vec_clusters.size(),0.0));
       for(size_t c=0; c< vec_clusters.size(); c++){
    //Loop over all hits in this clusters
    for(size_t c2 = 0; c2 < vec_clusters.size(); c2++){
    int n2_hits  = vec_clusters[c2].getHits().size(); 
    int n2_matched = 0; 
    if(vec_clusters[c2].getPlane() == vec_clusters[c].getPlane()) continue; //ignore clusters on same plane
    for(auto &h : vec_clusters[c].getHits()){
    //get time spread of this hit.
    double pp = h->PeakTimePlusRMS(1.0);
    double pm = h->PeakTimeMinusRMS(1.0);
    for(auto &h2 : vec_clusters[c2].getHits()){
    if(h2->PeakTime() < pp && h2->PeakTime() > pm) n2_matched++; 
    }
    }
    percent_matched[c][c2] = ((double)n2_matched)/((double)n2_hits)*100.0;
    std::cout<<" Cluster "<<c<<" on plane "<<vec_clusters[c].getPlane()<<" matches with "<<percent_matched[c][c2]<<" percent of hits of cluster "<<c2<<" on plane "<<vec_clusters[c2].getPlane()<<std::endl;
    }
    }
    */

    //This is where I will copy over a lot of the old SSS codebase
    std::vector<double> shr_start_3D= {vec_showers[0]->ShowerStart().X(), vec_showers[0]->ShowerStart().Y(),vec_showers[0]->ShowerStart().Z()};
    std::vector<double> shr_other_3D=  {vec_showers[0]->ShowerStart().X()+vec_showers[0]->Direction().X(),vec_showers[0]->ShowerStart().Y()+vec_showers[0]->Direction().Y(), vec_showers[0]->ShowerStart().Z()+vec_showers[0]->Direction().Z()};

    std::vector<std::vector<double>> shr_start_pt =   this->to2D(shr_start_3D);
    std::vector<std::vector<double>> shr_other_pt =   this->to2D(shr_other_3D);



    //Loop over all clusters
    for(size_t c=0; c< vec_clusters.size(); c++){

      auto hitz = vec_clusters[c].getHits(); // type of hitz: std::vector<art::Ptr<recob::Hit>>
      int num_hits_in_cluster = vec_clusters[c].getHits().size();
      int pl = vec_clusters[c].getPlane();

      //Need to modify this a bit
      auto ssscorz = SeaviewScoreCluster(pl,c+1, hitz ,vertex_chan[pl], vertex_tick[pl], vec_showers.front());
      vec_clusters[c].setScore(ssscorz);
      vec_clusters[c].setWireTickBasedLength(sqrt(pow((ssscorz.max_wire-ssscorz.min_wire)*wire_con, 2.0) + pow((ssscorz.max_tick - ssscorz.min_tick)*tick_con, 2.0)));
      BasicClusterCalorimetry(vec_clusters[c]);

      //This is just checking if its in, we can do this earlier; TODO
      //TODO: is_in_shower, get back to primary shower (at least available)
      //Sim Stuff
      //Draw Direction on plot
      //Delauney on here might be good, that said, we have a LOT of things. Hmm, cap at 50 hits maybe? 
      int is_in_shower = SeaviewCompareToShowers(pl,c+1, hitz ,vertex_chan[pl], vertex_tick[pl], vec_showers, showerToPFParticleMap, pfParticleToHitsMap,eps);

      std::string sname = "Cluster "+std::to_string(c)+" Hits: "+std::to_string(num_hits_in_cluster)+" PCA: "+std::to_string(ssscorz.pca_0)+" Theta: " +std::to_string(ssscorz.pca_theta)+" Wires: "+std::to_string(ssscorz.n_wires)+ " Ticks: "+std::to_string(ssscorz.n_ticks)+" ReMerged: "+std::to_string(is_in_shower);
      std::cout<<sname<<std::endl;

      if(ssscorz.pass){

        if(num_hits_in_cluster>0){


          TGraph * tmp = (TGraph*)vec_clusters[c].getGraph()->Clone(("tmp_"+std::to_string(pl)+std::to_string(c)).c_str());

          int Npts = 20;
          //TODO need to write this function
          TGraph * core  = (TGraph*)SeaviewGetNearestNpts(pl,c+1,hitz,vertex_chan[pl],vertex_tick[pl],Npts);

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

          //std::cout<<"Just Before Core Fit of "<<tmp->GetN()<<" pts between "<<chan_min[pl]<<" "<<chan_max[pl]<<" or "<<fmin<<" "<<fmax<<std::endl;

          double con;
          double slope;
          if(fmin==fmax){
            slope = 0;
            con = fmin;
          }else{
            core->Fit("pol1","Q","same",fmin,fmax); // fit to polynomial of degree 1, and plot it on the same pad
            con = core->GetFunction("pol1")->GetParameter(0);
            slope = core->GetFunction("pol1")->GetParameter(1);
          }

          double impact_parameter = 1e10;

          //rudimentary!
          for(double k=chan_min[pl]; k< chan_max[pl];k++){
            double y = slope*k+con;
            double dist = dist_point_point(k, y, vertex_chan[pl], vertex_tick[pl]);
            impact_parameter = std::min(impact_parameter,dist);
          }

          //Lets assume its potining back to the "vertex" and calculate a kinda_angle w.r.t to the shower
          //vertex_wire[i] vertex_tick[i] (already calcuated)
          //cluster closest point )ssscorz.close_wire and close_tick
          //recob::Shower start point, convered to wire tick.
          std::vector<double> vec_c = {(double)(vertex_chan[pl]-ssscorz.close_wire), (double)(vertex_tick[pl]-ssscorz.close_tick)};
          std::vector<double> vec_s = {(double)vertex_chan[pl]-shr_start_pt[pl][0], (double)vertex_tick[pl]-shr_start_pt[pl][1]};
          double l_c = sqrt(pow(wire_con*vec_c[0],2)+pow(vec_c[1]*tick_con,2));
          double l_s = sqrt(pow(wire_con*vec_s[0],2)+pow(vec_s[1]*tick_con,2));
          double kinda_angle = acos((wire_con*vec_s[0]*wire_con*vec_c[0]+vec_c[1]*vec_s[1]*tick_con*tick_con )/(l_c*l_s));


          std::cout<<"SSSNEW "<<this->tag<<std::endl;
          std::cout<<pl<<" Num Hits "<<num_hits_in_cluster<<std::endl;
          std::cout<<pl<<" Num Wires "<< (int)ssscorz.n_wires<<std::endl;
          std::cout<<pl<<" Num Ticks "<< (int)ssscorz.n_ticks<<std::endl;
          std::cout<<pl<<" Plane "<<pl<<std::endl;
          std::cout<<pl<<" PCA "<<ssscorz.pca_0<<std::endl;
          std::cout<<pl<<" Impact "<<impact_parameter<<std::endl;
          std::cout<<pl<<" Fit Slope "<<slope<<std::endl;
          std::cout<<pl<<" Fit Constant "<<con<<std::endl;
          std::cout<<pl<<" Mean Tick "<<ssscorz.mean_tick<<std::endl;
          std::cout<<pl<<" Max Tick "<<ssscorz.max_tick<<std::endl;
          std::cout<<pl<<" Min Tick "<<ssscorz.min_tick<<std::endl;
          std::cout<<pl<<" Mean Wire "<<ssscorz.mean_wire<<std::endl;
          std::cout<<pl<<" Max Wire "<<ssscorz.max_wire<<std::endl;
          std::cout<<pl<<" Min Wire "<<ssscorz.min_wire<<std::endl;
          std::cout<<pl<<" Kinda Angle "<<kinda_angle<<std::endl;

          vec_clusters[c].setShowerRemerge(is_in_shower);
          vec_clusters[c].setImpactParam(impact_parameter);
          vec_clusters[c].setFitSlope(slope);
          vec_clusters[c].setFitCons(con);
          vec_clusters[c].setAngleWRTShower(kinda_angle);
        }
      } // if ssscorz has passed

    }//cluster loop

    vec_in_clusters = vec_clusters;




    /*
    //Relative to showers!
    for(size_t s=0; s<vec_showers.size(); ++s){
    std::vector<double> shr_start_3D= {vec_showers[s]->ShowerStart().X(), vec_showers[s]->ShowerStart().Y(),vec_showers[s]->ShowerStart().Z()};
    std::vector<double> shr_other_3D=  {vec_showers[s]->ShowerStart().X()+vec_showers[s]->Direction().X(),vec_showers[s]->ShowerStart().Y()+vec_showers[s]->Direction().Y(), vec_showers[s]->ShowerStart().Z()+vec_showers[s]->Direction().Z()};

    std::vector<std::vector<double>> start_pt =   this->to2D(shr_start_3D);
    std::vector<std::vector<double>> other_pt =   this->to2D(shr_other_3D);

    for(int i=0; i<3; i++){
    double slope = (start_pt[i][1]-other_pt[i][1])/(start_pt[i][0]-other_pt[i][0]);
    double inter = start_pt[i][1]-slope*start_pt[i][0];

    double x1_plot = other_pt[i][0];//chan_min[i]-chan_shift;
    double y1_plot = slope*x1_plot+inter;

    double x2_plot;
    if(other_pt[i][0]<start_pt[i][0]){
    x2_plot = chan_max[i]+chan_shift;
    }else{
    x2_plot = chan_min[i]-chan_shift;    
    }
    double y2_plot = slope*x2_plot+inter;

    TLine *l = new TLine(x1_plot, y1_plot, x2_plot, y2_plot);


    //Ok now find the min distance from this line to all clusters on this plane

    for(size_t c=0; c< vec_clusters.size(); c++){
    int pl = vec_clusters[c].getPlane();
    if(pl != i ) continue;

    double min_d = 1e10;
    int min_p = -9;

    for(size_t p = 0; p< vec_clusters[c].getPts().size(); p++){
    std::vector<double> pt = (vec_clusters[c].getPts())[p];
    double dist = this->dist_line_point({x1_plot*wire_con,y1_plot*tick_con}, {x2_plot*wire_con,y2_plot*tick_con}, {pt[0]*wire_con, pt[1]*tick_con});

    if(dist< min_d){
    min_p = (int)p;
    min_d = dist;
    }
    }

    }



    }

    }
    */

    return {0};


  }

  double SEAviewer::dist_line_point( const std::vector<double>&X1, const std::vector<double>& X2, const std::vector<double>& point){
    // convert {wire, tick} coordinate to [cm, cm] coordinate
    double x1 =X1.at(0)*wire_con;
    double y1 =X1.at(1)*tick_con;

    double x2 =X2.at(0)*wire_con;
    double y2 =X2.at(1)*tick_con;

    double x0 =point.at(0)*wire_con;
    double y0 =point.at(1)*tick_con;

    double x10 = x1-x0;
    double y10 = y1-y0;

    double x21 = x2-x1;
    double y21 = y2-y1;

    double t = -(x10*x21+y10*y21)/fabs(x21*x21+y21*y21);

    double d2 = pow(x1-x0,2)+pow(y1-y0,2)+2*t*((x2-x1)*(x1-x0)+(y2-y1)*(y1-y0))+t*t*( pow(x2-x1,2)+pow(y2-y1,2));


    return sqrt(d2);

  }

  cluster_score SEAviewer::SeaviewScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower){
    cluster_score score(p,cl);
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
      double dd = dist_point_point(h_wire, h_tick, vertex_wire, vertex_tick);
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

  int SEAviewer::SeaviewCompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick, std::vector<art::Ptr<recob::Shower>>& showers, const std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showerToPFParticleMap, const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap, double eps){


    for(size_t s =0; s< showers.size(); s++){
      art::Ptr<recob::Shower> shower = showers[s];
      art::Ptr<recob::PFParticle> pfp = showerToPFParticleMap.at(shower);  //key has to be in the map, otherwise out-of-range error
      std::vector<art::Ptr<recob::Hit>> showerhits = pfParticleToHitsMap.at(pfp);

      bool in_primary_shower = false;
      for(size_t h = 0; h< hitz.size(); h++){
        auto hit = hitz[h];
        double h_wire = (double)hit->WireID().Wire;
        double h_tick = (double)hit->PeakTime();


        for(auto &sh: showerhits){

          if(sh->View() != hit->View()) continue;  //Guanqun: if not on the same plane?

          double sh_wire = (double)sh->WireID().Wire;
          double sh_tick = (double)sh->PeakTime();

          double dist = dist_point_point(sh_wire, sh_tick, h_wire, h_tick);

          if(dist<=eps){
            in_primary_shower = true;
            return (int)s;
          }

        }

      } // end of hitz loop

      if(in_primary_shower){
        return (int)s;
      }
    }


    return -1;
  }



  TGraph* SEAviewer::SeaviewGetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts){

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

      double dd = dist_point_point(h_wire, h_tick, vertex_wire, vertex_tick);
      all_wire.push_back(h_wire);   
      all_tick.push_back(h_tick);   
      all_dist.push_back(dd);
    }

    // sorted_in has indices of elements in all_dist in descending order
    std::vector<size_t> sorted_in = seaview_sort_indexes(all_dist);
    size_t max_e = std::min((size_t)Npts,hitz.size());

    for(size_t i =0; i<max_e; i++){ // get the position ({wire, tick}) of the closest 'max_e' points.
      t_wire.push_back(all_wire[sorted_in[hitz.size()-1-i]]);
      t_tick.push_back(all_tick[sorted_in[hitz.size()-1-i]]);
    }

    return new TGraph(t_wire.size(),&t_wire[0],&t_tick[0]);
  }


  void SEAviewer::BasicClusterCalorimetry(cluster& cl){

    //grab all hits in cluster
    const std::vector<art::Ptr<recob::Hit>>& hitz = cl.getHits();
    cl.f_ADC_hist.StatOverflows(kTRUE);

    for(auto& h : hitz){
      cl.f_ADC_hist.Fill(h->ROISummedADC());
    }

    cl.f_meanADC = cl.f_ADC_hist.GetMean();
    cl.f_ADC_RMS = cl.f_ADC_hist.GetRMS();
    return;
  }



  void SEAviewer::SetClusterLegend(int cluster, double energy, int is_matched, int matched_pdg, double overlay_fraction){

    //grab the plane number, and impact parameter of the cluster
    int plane = vec_clusters.at(cluster).getPlane();
    double min_ioc_to_shower  = vec_clusters.at(cluster).getMinHitIOC();

    //need to use stringstream to control the number of digits..
    std::ostringstream ss1, ss2, ss3;
    ss1 << std::setprecision(1) << std::fixed << energy;
    ss2 << std::setprecision(1) << std::fixed << min_ioc_to_shower;
    ss3 << std::setprecision(2) << std::fixed << overlay_fraction;

    std::string legend;
    //add the truth information to the legend if the cluster is matched to a MCParticle
    if(is_matched == 1){
      legend = "#splitline{" + std::to_string(plane) + ", " + ss1.str() + "MeV, Min IOC: " 
        + ss2.str() + "}{#splitline{Matched: " + (is_matched == 1 ? "true" : "false") +", PDG: " 
        + std::to_string(matched_pdg) + "}{Ovelay Frac: "+ ss3.str() + "}}";
    }
    else{
      legend = std::to_string(plane) + ", " + ss1.str() + "MeV, Min IOC: " + ss2.str();
    }
    vec_clusters.at(cluster).setLegend(legend);
  }


  void SEAviewer::format_legend(std::string &leg, double arg1, double arg2, double arg3){
    std::ostringstream ss1, ss2, ss3;
    ss1 << std::setprecision(1) << std::fixed << arg1;
    ss2 << std::setprecision(2) << std::fixed << arg2;
    ss3 << std::setprecision(1) << std::fixed << arg3;

    if(leg == "Shower"){
      leg = "#splitline{" + leg + ": " + ss1.str() + " MeV | " + ss2.str() + " cm }{conv. dist | "
        + ss3.str() + " impact par.}";
      //leg += ": " + ss1.str() + " MeV, " + ss2.str() + " cm conv. dist.";

    }else{
      //for tracks, 3rd argument is not used
      leg += ": "+ ss1.str() + " cm | " + ss2.str() + " PCA";
    }
  }
}
