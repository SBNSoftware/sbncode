#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "TGraph.h"

#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"

namespace single_photon
{

  struct sss_score{
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
    double close_wire; /* wire of hit that's closest to vertex */
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

    int n_wires; /* number of wires hits correspond to */
    int n_ticks;

    bool pass;

    sss_score(int ip, int cl): plane(ip), cluster_label(cl){};
  }; //end of class sss_score


    TGraph* GetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts);

    sss_score ScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower);

    int CompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,
            const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      const   std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap,                    double eps);



    std::vector<double>SecondShowerMatching(
      std::vector<art::Ptr<recob::Hit>>& hitz,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< int ,art::Ptr<simb::MCParticle>>  & MCParticleToTrackIdMap
  );


    //************************************************ Shower Search Slice Second SSS3D ********** /

    void SecondShowerSearch3D(
    std::vector<art::Ptr<recob::Shower>> & showers,
    std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & NormalShowerToPFParticleMap,  
    std::vector<art::Ptr<recob::Track>> & tracks, 
    std::map<art::Ptr<recob::Track>,  
    art::Ptr<recob::PFParticle>> & NormalTrackToPFParticleMap, 
    art::Event const & evt );

    void SimpleSecondShowerCluster();

    std::pair<bool, std::vector<double>> clusterCandidateOverlap(const std::vector<int> & candidate_indices, const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks);

   
  std::pair<int, std::pair<std::vector<std::vector<double>>, std::vector<double>>> GroupClusterCandidate(int num_clusters,  const std::vector<int>& cluster_planes, const std::vector<double>& cluster_max_ticks, const std::vector<double>& cluster_min_ticks);

//isolation.h
bool  map_max_fn(const std::pair<art::Ptr<recob::Hit>,double> p1, const std::pair<art::Ptr<recob::Hit>,  double> p2){
  return (p1.second < p2.second);
}

// override function of sorts for min_element function comparison
bool  map_min_fn(const std::pair<art::Ptr<recob::Hit>,double> p1, const std::pair<art::Ptr<recob::Hit>,  double> p2){
  return (p1.second > p2.second);
}

  void IsolationStudy(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Track>>& tracks, 
      const std::vector<art::Ptr<recob::Shower>>& showers, 
      detinfo::DetectorPropertiesData const & theDetector);


}
