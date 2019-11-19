#ifndef __sbnanalysis_CURS_HH
#define __sbnanalysis_CURS_HH

#include <array>

#include "TVector3.h"

#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "../Data/RecoEvent.h"

namespace ana {
 namespace SBNOsc {


class Cuts {
public:
  static const unsigned nCuts = 9; //!< total number of cuts
  static const unsigned nTruthCuts = 5; //!< Total number of truth cuts
  /**
 * Initialize this class.
 * \param cfg fhicl configuration for the class
 * \geometry A pointer to the geoemtry of the detector to configure
 */
  void Initialize(const fhicl::ParameterSet &cfg, const geo::GeometryCore *geometry);
  /** 
 * Process each cut associated with reconstructed events
 * \param event The reconstructed event information
 * \param reco_vertex_index The index of the candidate reconstructed neutrino vertex into the list of such vertices in the RecoEvent
 *
 * \return A list of bool's of whether the reco event passes each cut
 */
  std::array<bool, nCuts> ProcessRecoCuts(const numu::RecoEvent &event, unsigned reco_vertex_index) const;

  /**
 * Process each cut associated with true events
 * \param event Event information
 * \param truth_vertex_index The index of the true neutrino vertex into the lsit of true interactions in the RecoEvent
 * 
 * \return A list of bool's of whether the true event passes each cut
 */
  std::array<bool, nTruthCuts> ProcessTruthCuts(const numu::RecoEvent &event, unsigned truth_vertex_index) const;

  /**
 * Select a reco event based on the cut values provided by ProcessRecoCuts
 * \param cuts the list of cuts returned by ProcessRecoCuts
 *
 * \return whether to select this reconstructed neutrino vertex candidate
 */
  bool SelectReco(std::array<bool, nCuts> &cuts);
  /**
 * Test whether a point is in the configured fiducial volume
 * \param v The point to test
 *
 * \return Whether the point is in the configured fiducial volume
 */
  bool InFV(const TVector3 &v) const;
  /**
 * Test whether a point is in the configured fiducial volume
 * \param v The point to test
 *
 * \return Whether the point is in the configured fiducial volume
 */
  bool InFV(const geo::Point_t &v) const;
  /**
 * Test whether a point is in the configured track containment volume
 * \param v The point to test
 *
 * \return Whether the point is in the configured track containment volume
 */
  bool InCosmicContainment(const TVector3 &v) const;

  /**
 * Test whether a point is in the configured track containment volume
 * \param v The point to test
 *
 * \return Whether the point is in the configured track containment volume
 */
  bool InCalorimetricContainment(const TVector3 &v) const;

  /**
 * Gets the time of the CRT match to a track.
 * \param track The track object
 *
 * \return The time [us] of the CRT match. Returns nonsense if no such match exists.
 */
  float CRTMatchTime(const numu::RecoTrack &track) const;

  /**
 * Whether a TPC track has a CRT hit match
 * \param track the TPC track
 * \return Whether a TPC track has a CRT hit match
 */
  bool HasCRTHitMatch(const numu::RecoTrack &track) const;

  /**
 * Returns whether a time value is within the configured beam spill window
 * \param time The time to test
 * \return True when the time is in the beam spill window.
 */
  bool TimeInSpill(float time) const;

  /**
 * Whether a TPC track has a CRT track match
 * \param track the TPC track
 * \return Whether a TPC track has a CRT track match
 */
  bool HasCRTTrackMatch(const numu::RecoTrack &track) const;
  

private:
  struct VolYZ {
    std::array<double, 2> Y;
    std::array<double, 2> Z;
  };

  struct Config {
    double trackMatchCompletionCut;
    float TruthCompletion;
    float TruthMatchDist;
    float CRTHitDist;
    std::array<float, 2> CRTHitTimeRange;
    float CRTTrackAngle;
    float TrackLength;
    float MCSTrackLength;
    std::vector<geo::BoxBoundedGeo> fiducial_volumes;
    std::vector<VolYZ> cosmic_containment_volumes;
    std::vector<geo::BoxBoundedGeo> active_volumes;
    std::vector<geo::BoxBoundedGeo> calorimetric_containment_volumes;
    bool TruthFlashMatch;
  };

  Config fConfig;

};

  }
}

#endif

