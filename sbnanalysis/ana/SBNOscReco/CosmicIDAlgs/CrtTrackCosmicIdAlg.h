#ifndef CRTTRACKCOSMICIDALG_H_SEEN
#define CRTTRACKCOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// CrtTrackCosmicIdAlg.h
//
// Functions for CRTTrack match cosmic tagger
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// sbndcode
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"

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

  class CrtTrackCosmicIdAlg {
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

      fhicl::Table<sbnd::CRTTrackMatchAlg::Config> TrackMatchAlg {
        Name("TrackMatchAlg"),
        Comment("")
      };

      fhicl::Table<BeamTime> BeamTimeLimits {
        Name("BeamTimeLimits"),
        Comment("")
      };

    };

    CrtTrackCosmicIdAlg(const Config& config);

    CrtTrackCosmicIdAlg(const fhicl::ParameterSet& pset) :
      CrtTrackCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

    CrtTrackCosmicIdAlg();

    ~CrtTrackCosmicIdAlg();

    void reconfigure(const Config& config);

    // Tags track as cosmic if it matches a CRTTrack
    bool CrtTrackCosmicId(recob::Track track, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event);

    // Getter for matching algorithm
    sbnd::CRTTrackMatchAlg TrackAlg() const {return trackMatchAlg;}

  private:

    sbnd::CRTTrackMatchAlg trackMatchAlg;
    double fBeamTimeMin;
    double fBeamTimeMax;

  };

}

#endif
