////////////////////////////////////////////////////////////////////////
// Class:       AnodeCrossT0
// Plugin Type: producer (Unknown Unknown)
// File:        AnodeCrossT0_module.cc
//
// Generated at Sun May 29 12:26:30 2022 by Gray Putnam using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <memory>

class AnodeCrossT0;


class AnodeCrossT0 : public art::EDProducer {
public:
  explicit AnodeCrossT0(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnodeCrossT0(AnodeCrossT0 const&) = delete;
  AnodeCrossT0(AnodeCrossT0&&) = delete;
  AnodeCrossT0& operator=(AnodeCrossT0 const&) = delete;
  AnodeCrossT0& operator=(AnodeCrossT0&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  art::InputTag fPFPproducer;
  art::InputTag fTRKproducer;
  double fMinTimeTickInset;
  double fMaxTimeTickInset;
  double fInsetY;
  double fInsetZ;
  double fLengthCut;
  bool fSilenceMissingDataProducts;
};


AnodeCrossT0::AnodeCrossT0(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fPFPproducer(p.get<art::InputTag>("PFPproducer")),
    fTRKproducer(p.get<art::InputTag>("TRKproducer")),
    fMinTimeTickInset(p.get<double>("MinTimeTickInset")),
    fMaxTimeTickInset(p.get<double>("MaxTimeTickInset")),
    fInsetY(p.get<double>("InsetY")),
    fInsetZ(p.get<double>("InsetZ")),
    fLengthCut(p.get<double>("LengthCut")),
    fSilenceMissingDataProducts(p.get<bool>("SilenceMissingDataProducts", false))
{
  produces<std::vector<anab::T0>>();
  produces<art::Assns<recob::PFParticle, anab::T0>>();

}

void AnodeCrossT0::produce(art::Event& e)
{
  // output
  std::unique_ptr<std::vector<anab::T0>> outT0(new std::vector<anab::T0>);
  std::unique_ptr<art::Assns<recob::PFParticle, anab::T0>> assn(new art::Assns<recob::PFParticle, anab::T0>);

  art::PtrMaker<anab::T0> t0PtrMaker {e};

  // Services
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const dprop =
    art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  double maxy = geometry->TPC(0, 0).ActiveBoundingBox().MaxY();
  double maxz = geometry->TPC(1, 0).ActiveBoundingBox().MaxZ();
  double minz = geometry->TPC(0, 0).ActiveBoundingBox().MinZ();

  double maxtick = dprop.NumberTimeSamples() - fMaxTimeTickInset;

  // Reconstructed Information
  std::vector<art::Ptr<recob::PFParticle>> PFParticleList;
  try {
    art::ValidHandle<std::vector<recob::PFParticle>> pfparticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::fill_ptr_vector(PFParticleList, pfparticles);
  }
  catch(...) {
      std::cout << "PFP's with tag: " << fPFPproducer << " not present.\n";
      // Real data may have missing products -- just ignore the event
      if (fSilenceMissingDataProducts) return;
      else throw;
  }

  // T0 -- associated data
  art::FindManyP<anab::T0> fmT0(PFParticleList, e, fPFPproducer);

  // Track - associated data
  art::FindManyP<recob::Track> fmTracks(PFParticleList, e, fTRKproducer);

  for (unsigned i_pfp = 0; i_pfp < PFParticleList.size(); i_pfp++) {
    // Only look for cosmic-id'd tracks 
    if (!PFParticleList[i_pfp]->IsPrimary()) continue;

    // Ignore ones that crossed the Cathode and already got a T0
    std::vector<art::Ptr<anab::T0>> maybe_t0 = fmT0.at(i_pfp);
    if (!maybe_t0.empty()) continue;

    std::vector<art::Ptr<recob::Track>> maybe_track = fmTracks.at(i_pfp);
    if (maybe_track.empty()) continue;

    const recob::Track &trk = *maybe_track[0];
    
    geo::Point_t start = trk.Start();
    geo::Vector_t start_dir = trk.StartDirection();

    bool anode_crosser = start_dir.Y() < 0 && (start.Y() < maxy - fInsetY) \
        && (start.Z() < maxz - fInsetZ) && (start.Z() > minz + fInsetZ) \
        && trk.Length() > fLengthCut;

    if (!anode_crosser) continue;

    // Get the hit information
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmtrkHits(maybe_track, e, fTRKproducer);
    const std::vector<art::Ptr<recob::Hit>> &hits = fmtrkHits.at(0);
    const std::vector<const recob::TrackHitMeta *> &thms = fmtrkHits.data(0);

    // Get the first and last time hit on the track (use the collection plane only)
    double min_peak;
    int i_hit0 = -1;
    for (unsigned i_thm = 0; i_thm < thms.size(); i_thm++) {
      const recob::Hit &hit = *hits[i_thm];
      bool badhit = (thms[i_thm]->Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!trk.HasValidPoint(thms[i_thm]->Index()));

      if (!badhit && hit.SignalType() == geo::kCollection) {
        double peaktime = hit.PeakTime();
        if (i_hit0 < 0 || peaktime < min_peak) {
          i_hit0 = (int) i_thm;
          min_peak = peaktime;
        }
      }
    }

    double max_peak;
    int i_hitl = -1;
    for (unsigned i_thm = 0; i_thm < thms.size(); i_thm++) {
      const recob::Hit &hit = *hits[i_thm];
      bool badhit = (thms[i_thm]->Index() == std::numeric_limits<unsigned int>::max()) ||
                    (!trk.HasValidPoint(thms[i_thm]->Index()));

      if (!badhit && hit.SignalType() == geo::kCollection) {
        double peaktime = hit.PeakTime();
        if (i_hitl < 0 || peaktime > max_peak) {
          i_hitl = (int) i_thm;
          max_peak = peaktime;
        }
      }
    }

    if (i_hit0 > 0 && i_hitl > 0) {
      const recob::Hit &hit0 = *hits[i_hit0];
      geo::Point_t p0 = trk.LocationAtPoint(thms[i_hit0]->Index());
      double peak0 = hit0.PeakTime();

      const recob::Hit &hitl = *hits[i_hitl];
      geo::Point_t pl = trk.LocationAtPoint(thms[i_hitl]->Index());
      double peakl = hitl.PeakTime();

      // Do checks on the hits
      bool anode_crosser_hit = (peak0 > fMinTimeTickInset) && (peakl < maxtick) && (p0.Y() > pl.Y());

      if (anode_crosser_hit) {
        // All done! Get the t0
        double t0 = clock_data.TPCTick2TrigTime(peak0)*1e3; // us -> ns

        anab::T0 aT0(t0, 0, 0);
        outT0->push_back(aT0);
        art::Ptr<anab::T0> thisT0Ptr = t0PtrMaker(outT0->size() - 1);
        assn->addSingle(PFParticleList[i_pfp], thisT0Ptr);
      }
    }
  }

  e.put(std::move(outT0));
  e.put(std::move(assn));

}

DEFINE_ART_MODULE(AnodeCrossT0)
