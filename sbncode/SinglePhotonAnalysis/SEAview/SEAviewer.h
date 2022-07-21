#ifndef SEAVIEWER_H
#define SEAVIEWER_H

#include <iostream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include <string>

#include <memory>

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPrincipal.h"
#include "TFitResultPtr.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TH1D.h"
#include "TEllipse.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>


#include "seaDBSCAN.h"
namespace seaview {

  class SEAviewer;

  template <typename T>
    std::vector<size_t> seaview_sort_indexes(const std::vector<T> &v) {

      std::vector<size_t> idx(v.size());
      std::iota(idx.begin(), idx.end(), 0); //fill the range with sequentially increasing values

      // sort indexes based on comparing values in v (descending order)
      std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

      return idx;
    }



  struct cluster_score{
    int plane;
    int cluster_label;
    double point_score;
    int n_hits;

    double mean_wire;
    double max_wire;
    double min_wire;
    double mean_tick;
    double max_tick;
    double min_tick;

    double close_tick;
    double close_wire;
    double angle;//w.r.t shower primary

    double impact_parameter;

    double max_dist_tick;
    double mean_dist_tick;
    double min_dist_tick;
    double max_dist_wire;
    double mean_dist_wire;
    double min_dist_wire;

    double mean_dist;
    double max_dist;
    double min_dist;

    double pca_0;
    double pca_1;
    double pca_theta;

    int n_wires; /* number of wires spanned by the cluster */
    int n_ticks;  /* number of ticks spanned by the cluster */

    bool pass;

    cluster_score(int ip, int cl): plane(ip), cluster_label(cl){};
  };  //end of cluster_score class




  class cluster {
    friend class SEAviewer;

    public:

    cluster(int ID, int plane, std::vector<std::vector<double>> &pts, std::vector<art::Ptr<recob::Hit>> &hits) :f_ID(ID), f_plane(plane), f_pts(pts), f_hits(hits), f_score(0,0), f_shower_remerge(-1){

      f_npts = f_pts.size();
      if(pts.size() != hits.size()){
        std::cerr<<"seaviewer::cluster, input hits and pts not alligned"<<std::endl;
      }
      std::vector<double> wires(f_npts);
      std::vector<double> ticks(f_npts);
      for(int p =0; p< f_npts; ++p){
        wires[p]=f_pts[p][0];
        ticks[p]=f_pts[p][1];
      }
      TGraph af_graph(f_npts, &wires[0], &ticks[0]);
      f_graph = af_graph;

    };
    int setScore(cluster_score &in_score){ f_score = in_score;return 0;}
    void setLegend(const std::string &in_leg){
      f_legend = in_leg;
    }

    void setWireTickBasedLength( double d) { f_wire_tick_based_length = d;}
    double getWireTickBasedLength() const {return f_wire_tick_based_length; }
    bool isTrackAnalyzed() const {return f_track_treated; }
    cluster_score * getScore(){return &f_score;};
    int getID() const {return f_ID;}
    int getPlane() const {return f_plane; }
    std::vector<std::vector<double>> getPTS() const {return f_pts;}
    TGraph * getGraph(){ return &f_graph;}
    TH1D* getADCHist() {return &f_ADC_hist;}
    const TGraph * getGraph() const { return &f_graph;}
    const std::string &getLegend() const {return f_legend; }
    const std::vector<art::Ptr<recob::Hit>>&  getHits(){return f_hits;}
    int getShowerRemerge() const {return f_shower_remerge;}
    int setShowerRemerge(int remerge_in){
      f_shower_remerge = remerge_in;
      return f_shower_remerge;
    }
    double getMeanADC() const { return f_meanADC; }
    double getADCrms() const {return f_ADC_RMS;}

    //----------------------- second-shower relatd function -----------------------
    void setImpactParam(double d) {f_ImpactParameter = d; } 

    /* brief: set the angle between direction of second shower candidate cluster and direction of the primary shower*/
    void setAngleWRTShower(double d) {f_AngleWRTShower = d;}
    void setFitSlope(double d) { f_FitSlope = d;}
    void setFitCons(double d) {  f_FitCons = d;}
    double getAngleWRTShower() const {return f_AngleWRTShower;}
    double getFitSlope() const {return f_FitSlope; }
    double getFitCons() const {return f_FitCons;}
    double getImpactParam() const {return f_ImpactParameter; } 


    // ----------------------- track search related function -----------------------
    double getMinHitImpactParam() const {return f_min_impact_parameter_to_shower; }
    double getMinHitConvDist() const { return f_min_conversion_dist_to_shower_start; }
    double getMinHitIOC() const {return f_min_ioc_to_shower_start;}
    double getIOCbasedLength() const {return f_ioc_based_length; }
    size_t getTrackStartIdx() const {return start_hit_idx; }
    size_t getTrackEndIdx() const {return end_hit_idx;}
    double getMeanADCFirstHalf() const { return f_mean_ADC_first_half; }
    double getMeanADCSecondHalf() const {return f_mean_ADC_second_half; }
    double getMeanADCRatio() const {return f_mean_ADC_first_to_second_ratio; }
    double getTrackAngleToShowerDirection() const {return f_angle_wrt_shower_direction; }
    double getLinearChi() const {return f_fit_chi2; }

    /* brief: check if this cluster is fully in given slice or not
     * return: 1 -> Fully in given Slice;   -1 --> Fully not in given slice;  0: parts in given slices, parts not
     */
    int InNuSlice(const std::map<int, std::vector<art::Ptr<recob::Hit>> >& sliceIDToHitsMap, int nuSliceID);

    // determine if the cluster is within the plot range
    // tick_max, tick_min, wire_max, and wire_min are the edges of the X axis(wire) and Y axis(tick)
    bool InRange(double tick_max, double tick_min, double wire_max, double wire_min) const{
      return f_score.min_wire < wire_max && f_score.max_wire > wire_min && f_score.max_tick > tick_min && f_score.min_tick < tick_max;
    }


    private:
    int f_ID;
    int f_npts;
    int f_plane;
    bool f_track_treated = false;  //boolean indicating whether the hits have been analyzed as track candidate

    std::vector<std::vector<double>> f_pts; //vector of {wire, tick} pairs of all the hits
    std::vector<art::Ptr<recob::Hit>> f_hits;
    std::vector<int> f_hit_group;   //group the hits in two groups: first half, second half (directionality-wise)
    cluster_score f_score;
    int f_shower_remerge = -1;  //index of the reco shower if the cluseter is close enough to a reco shower, otherwise -1.
    TGraph f_graph;        //2D {wire, tick} graph
    TH1D f_ADC_hist;    // histograms of ADC of every hit
    std::string f_legend; //legend of the f_graph


    // -------- track-like properties -------------------------

    //add a few parameters that are useful to find tracks in no recob::track events
    double f_min_impact_parameter_to_shower = 1e10;  // mininum impact parameter of hits to the recob::shower direction
    // will be default value if the cluster didn't pass cut on ssscore
    double f_min_conversion_dist_to_shower_start = 1e10;   //minimum distance of hits to the recob::shower start
    double f_min_ioc_to_shower_start = 1e10;               //minimum ioc of all hits to the recob::shower direction
    double f_ioc_based_length = -1.0;  // length of the cluster, calculated based on the IOC of hits
    double f_wire_tick_based_length = -1.0;

    size_t start_hit_idx = SIZE_MAX; //index of the start hit
    size_t end_hit_idx = SIZE_MAX;   //index of the end hit
    double f_mean_ADC_first_half = 0.0;
    double f_mean_ADC_second_half = 0.0;
    double f_mean_ADC_first_to_second_ratio = 0.0;

    double f_angle_wrt_shower_direction = -1.0; // angle between the cluster direction and the shower direction, in radian
    // for track search, for proton track veto
    double f_fit_chi2 = 0.0; //chi2 of the linear fit to the cluster (2D) 
    double f_ADC_RMS = -1.0;  //RMS of the summed ADC of hits
    double f_meanADC = -1.0;  // mean of hits ADC


    //-------- shower-like properties -------------------------------

    double f_ImpactParameter = -1.0; //impact parameter of the vertex wrt to the cluster
    double f_FitSlope = 0.0; //slope of the fitted shower/cluster direction
    double f_FitCons = 0.0;  //intercept of the fitted shower/cluster direction

    double f_AngleWRTShower = -1.0; //angle between cluster-vertex direction and primary_shower_start-vertex direction, assuming cluster and primary shower both point back to the vertex
    // specific for second shower search for 1g1p analysis

  };  // end of class cluster


  class SEAviewer {

    public:

      //constructor
      //Keng, DetectorProperties--> DetectorPropertiesData 
      SEAviewer(std::string tag,geo::GeometryCore const * geom,detinfo::DetectorPropertiesData const & theDetector );

      void configure(const fhicl::ParameterSet& pset){};

      int loadVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z);
      int addTrueVertex(double x, double y,double z);

      /* @brief: add all the given hits to be considered for clustering 
       * @note: these hits are subject to further selection by filterConsideredHits() before clustering
       */
      int addHitsToConsider(std::vector<art::Ptr<recob::Hit>>& hits);
      int addAllHits(std::vector<art::Ptr<recob::Hit>>& hits);

      /* @brief: filter given hits before further clustering - only use hits near vertex for clustering
       * @param: dist_to_vertex - distance to vertex in cm
       */
      int filterConsideredHits(double dist_to_vertex);    

      /* @brief: add hits for PFParticles
       * @param: leg - legend for this PFParticle, for plotting purposes
       * @param: if leg is 'shower', arg1 expects shower energy; arg2 expects conversion distance; arg3 expects impact parameter
       * @param: if leg is 'track', arg1 expects track length; arg2 expects PCA; arg3 currently not used
       */
      int addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits, std::string leg , double arg1 = 0.0, double arg2 = 0.0, double arg3=0.0);
      int setBadChannelList(std::vector<std::pair<int,int>> &in);
      int addShower(art::Ptr<recob::Shower>&shr);
      int addTrack(art::Ptr<recob::Track>&trk);
      std::vector<int> calcUnassociatedHits();
      int setHitThreshold(double);
      int Print(double plot_distance);
      int runseaDBSCAN(double min_pts, double eps);

      double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
        //WireCoordinate returns the index of the nearest wire to the specified position.
        double wire = geo.WireCoordinate(Y, Z, plane, fTPC, fCryostat);
        return wire;
      }

      double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorPropertiesData const& detprop){
        double time = detprop.ConvertXToTicks(X, plane, fTPC,fCryostat);
        return time;
      }

      /* @brief: given a 3D space point, calculate the [wire, tick] of the point on 3 planes */
      std::vector<std::vector<double>> to2D(std::vector<double> & threeD);


      /* @brief: given two points on a line, and another point (in 2D space), calculate the impact parameter 
       * @param: 3 parameter of std::vector<double> are {wire, tick} vectors 
       */
      double dist_line_point( const std::vector<double>&X1, const std::vector<double>& X2, const std::vector<double>& point);


      //distance between two {wire, tick} points, in cm
      double dist_point_point(double w1, double t1, double w2, double t2) const{
        return  sqrt(pow(w1*wire_con-w2*wire_con,2)+pow(t1*tick_con-t2*tick_con,2));
      }


      //@brief: given a cluster, loop through all its hits, calc and update info on mean ADC, ADC RMS.
      void BasicClusterCalorimetry(cluster& cl);

      // @brief: given primary shower, analyze all the clusters in a track-like way with the assumtion that primary shower will point back to the cluster
      // @param: vec_c: vector of clusters to be filled 
      std::vector<double> analyzeTrackLikeClusters(double eps, const std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap, std::vector<seaview::cluster> &vec_c );


      // @brief: given primary shower, analyze all the clusters, draw them on the canvas, together with fitted direction of the cluseter (with the assumption that the cluster and primary shower both pointing back to the vertex)
      // @param: vec_c: vector of clusters to be filled 
      std::vector<double> analyzeShowerLikeClusters(double eps, const std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap, std::vector<seaview::cluster> &vec_c );


      // @brief: analyze cluster cl as a track candidate, and save track-related information in the cluster object
      // @param: shower_start_pt_2D, shower_other_pt_2D: {wire, tick} coordinate of shower start, and another point on shower direction line, projected to the plane cl is on.
      void TrackLikeClusterAnalyzer(cluster &cl, const std::vector<double> &shower_start_pt_2D, const std::vector<double> &shower_other_pt_2D);


      // @brief: check if there is a hit in hitz close enought to one of the reco showers, if so return the index of that reco shower
      // @param: hitz is usually a cluster of unassociated hits
      int SeaviewCompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,   std::vector<art::Ptr<recob::Shower>>& showers, const std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showerToPFParticleMap, const std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap, double eps);

      // Guanqun: @brief: analyze a cluster of hits and summarize its property into an cluster_score 
      // argument shower is not used in the function body
      cluster_score SeaviewScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower);

      //@ brief: as name says, get Npts hits that're nearest from vertex
      // return the {wire, tick} info of these hits as a TGraph
      TGraph* SeaviewGetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts);

      void SetClusterLegend(int cluster, double energy, int is_matched, int matched_pdg, double overlay_fraction);



      //conversion from wire, tick to cm
      static constexpr double wire_con = 0.3;
      static constexpr double tick_con = 1.0/25.0;

    protected:

      int n_pfps;    // num of PFParticles (including shower & track)
      int n_showers; //num of showers
      int n_tracks;  //num of tracks

      std::string tag;
      double hit_threshold;
      bool has_been_clustered;
      // PFP, Plane: index  
      std::vector<std::vector<TGraph>> vec_graphs; //vector of graphs of [tick vs. wire] for hits of PFParticles.

      std::vector<std::string> vec_pfp_legend; //legend for each PFParticle, for plotting purpose
      // PFP, Plane: index
      std::vector<std::vector<std::vector<double>>> vec_ticks; //vector of ticks on each plane for all PFParticles
      std::vector<std::vector<std::vector<double>>> vec_chans; //vector of wires on each plane for all PFParticle

      geo::GeometryCore const * geom;
      detinfo::DetectorPropertiesData const & theDetector ;

      double tick_shift;
      double chan_shift;

      double tick_max; //min, max tick of all hits
      double tick_min;
      std::vector<double> chan_max; //min, max wire of all (including vertex)
      std::vector<double> chan_min;

      std::vector<std::pair<int,int>> m_bad_channel_list;

      //Vertex, size of 3 (on 3 planes)
      std::vector<double> vertex_tick; 
      std::vector<double> vertex_chan; 
      std::vector<TGraph> vertex_graph;

      bool plot_true_vertex;
      //True vertex, size of 3
      std::vector<double> true_vertex_tick; 
      std::vector<double> true_vertex_chan; 
      std::vector<TGraph> true_vertex_graph;

      //std::vector<art::Ptr<recob::Hit>> considered_hits; //all hits considered for clustering
      //std::vector<art::Ptr<recob::Hit>> all_hits;
      std::map<art::Ptr<recob::Hit>,bool> map_unassociated_hits;
      std::map<art::Ptr<recob::Hit>, bool> map_considered_hits;

      //Plane: index
      std::vector<TGraph> vec_unass_graphs; //graph of [tick vs wire] for unassociated hits that pass the hit threshold
      std::vector<std::vector<double>> vec_unass_ticks; //tick of unassso hits that pass threshold
      std::vector<std::vector<double>> vec_unass_chans;
      std::vector<std::vector<std::vector<double>>> vec_unass_pts; // [wire, tick] pair for unassociatd hits that pass threshold on each plane
      std::vector<std::vector<art::Ptr<recob::Hit>>> vec_unass_hits; //vector of unasso hits that pss hit threshold on each plane


      //Plane: index
      std::vector<TGraph> vec_all_graphs;  //graph of [tick vs wire] for all hits that are not in the slice
      std::vector<std::vector<double>> vec_all_ticks;  //tick of all hits that are not in the slice (grouped by plane #)
      std::vector<std::vector<double>> vec_all_chans;  //wire of all hits that are not in the slice on each plane.

      std::vector<int> num_clusters; //number of clusters for unassociated hits on each plane
      std::vector<std::vector<int>> cluster_labels; //one-to-one mapped cluster labels for unassociated hits in `vec_unass_pts`
      TRandom3 *rangen;

      // all clusters on all 3 planes, each cluster includes points/hits identified for that cluster
      std::vector<seaview::cluster> vec_clusters;  
      std::vector<art::Ptr<recob::Shower>> vec_showers; //vector of recob::Shower contained in this class
      std::vector<art::Ptr<recob::Track>> vec_tracks;

      //-----helper function-----------

      // form legend for recob::shower and recob::track objects
      void format_legend(std::string &leg, double arg1 = 0.0, double arg2 = 0.0, double arg3 = 0.0);

  };

  //define wire conversion, tick conversion factor
  constexpr double SEAviewer::wire_con;
  constexpr double SEAviewer::tick_con;

}// namespace

#endif

