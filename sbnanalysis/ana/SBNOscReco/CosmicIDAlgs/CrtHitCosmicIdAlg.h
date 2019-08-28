#ifndef CRTHITCOSMICIDALG_H_SEEN
#define CRTHITCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// CrtHitCosmicIdAlg.h
//
// Functions for CRTHit match cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// sbndcode
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"

// c++
#include <vector>


namespace ana{

  class CrtHitCosmicIdAlg {
  public:
  
    struct BeamTime {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> BeamTimeMin {
        Name("BeamTimeMin"),
        Comment("Minimum t0 tagged time cut")
      };

      fhicl::Atom<double> BeamTimeMax {
        Name("BeamTimeMax"),
        Comment("Maximum t0 tagged time cut")
      };

    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Table<sbnd::CRTT0MatchAlg::Config> T0Alg {
        Name("T0Alg"),
        Comment("Configuration for CRTHit matching")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("Configuration for t0 cut limits")
      };

    };

    CrtHitCosmicIdAlg(const Config& config);

    CrtHitCosmicIdAlg(const fhicl::ParameterSet& pset) :
      CrtHitCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    CrtHitCosmicIdAlg();

    ~CrtHitCosmicIdAlg();

    void reconfigure(const Config& config);

    // Returns true if matched to CRTHit outside beam time
    bool CrtHitCosmicId(recob::Track track, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event);

    // Getter for matching algorithm
    sbnd::CRTT0MatchAlg T0Alg() const {return t0Alg;}

  private:

    sbnd::CRTT0MatchAlg t0Alg;
    double fBeamTimeMin;
    double fBeamTimeMax;

  };

}

#endif
