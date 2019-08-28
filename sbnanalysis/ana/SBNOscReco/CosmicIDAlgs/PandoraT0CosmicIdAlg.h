#ifndef PANDORAT0COSMICIDALG_H_SEEN
#define PANDORAT0COSMICIDALG_H_SEEN


///////////////////////////////////////////////
// PandoraT0CosmicIdAlg.h
//
// Functions for pandora t0 cosmic tagger
// Pandora removed these particles anyway
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"

// c++
#include <vector>


namespace ana{

  class PandoraT0CosmicIdAlg {
  public:

    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("")
      };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> PandoraLabel {
        Name("PandoraLabel"),
        Comment("")
      };

      fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
        Name("TpcTrackModuleLabel"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    PandoraT0CosmicIdAlg(const Config& config);

    PandoraT0CosmicIdAlg(const fhicl::ParameterSet& pset) :
      PandoraT0CosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    PandoraT0CosmicIdAlg();

    ~PandoraT0CosmicIdAlg();

    void reconfigure(const Config& config);

    // Finds any t0s associated with track by pandora, tags if outside beam
    bool PandoraT0CosmicId(recob::Track track, const art::Event& event);

    // Finds any t0s associated with pfparticle by pandora, tags if outside beam
    bool PandoraT0CosmicId(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event);

  private:

    art::InputTag fPandoraLabel;
    art::InputTag fTpcTrackModuleLabel;
    double fBeamTimeMin;
    double fBeamTimeMax;

  };

}

#endif
