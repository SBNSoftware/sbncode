#ifndef SBNCODE_SINGLEPHOTONANALYSIS_LIBARIES_PROCESSORS_H
#define SBNCODE_SINGLEPHOTONANALYSIS_LIBARIES_PROCESSORS_H

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace single_photon
{
    int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input);

    int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected);

    int spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected);

    void CollectMCParticles(
    const art::Event &evt,
    const std::string &label,
    std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,
    std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> &particlesToTruth,
    std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap);

    void CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector);

  void BuildMCParticleHitMaps(
      const art::Event &evt, 
      const std::string &label, 
      const std::vector<art::Ptr<recob::Hit>> &hitVector,   
      std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         
      std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, 
      const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, 
      std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap);

    bool Pi0PreselectionFilter();

    bool Pi0PreselectionFilter2g0p();

    bool IsEventInList(int run, int subrun, int event);

//----------- Above are migrated from Singlephoton_module.cc

    //determines if a point is inside the rectangle by summing the areas of the four triangles made by 
    //if the point is inside, the sum of the triangles should exactly equal the area of the rectangle
    //also returns true if the point is on the boundary
    bool isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);


  //helpers for calculating calometry
    double CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits);

    double CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int this_plane);

    double GetQHit(art::Ptr<recob::Hit> thishitptr, int plane);


    double QtoEConversion(double Q);


    std::vector<double> CalcdEdxFromdQdx(std::vector<double> dqdx);


    std::vector<double> CalcdQdxShower(
            const art::Ptr<recob::Shower>& shower,
            const std::vector<art::Ptr<recob::Cluster>> & clusters, 
            std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane,
      double triggeroffset,
      detinfo::DetectorPropertiesData const & theDetector);

  std::vector<double> getPitch(TVector3 shower_dir);


    double getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);



    int getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);



    double triangle_area(double a1, double a2, double b1, double b2, double c1, double c2);

    int quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area);

    int delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area);

}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_LIBARIES_PROCESSORS_H
