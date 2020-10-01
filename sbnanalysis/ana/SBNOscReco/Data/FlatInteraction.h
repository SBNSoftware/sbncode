#ifndef _sbnumurecodata_FlatInteraction_hh
#define _sbnumurecodata_FlatInteraction_hh

namespace numu::flat {

struct TrackTruth {
  float momentum[3];
  float wall_enter;
  float wall_exit;
  float time;
  float is_contained;
  float pdgid;
  float is_cosmic;
};

struct PrimaryTrack {
  float length;
  float costh;
  float range_momentum;
  float mcs_momentum;
  float crt_hit_distance;
  float crt_hit_time;
  float crt_track_angle;
  float crt_track_time;
  float start[3];
  float end[3];
  TrackTruth truth;
};

struct TrueNeutrino {
  float E;
  float Q2;
  float vertex[3];
};

struct EventInfo {
  float crt_hit_times[10];
  float crt_hit_pes[10];
  float pass_trig;
};

struct EventMeta {
  float n_gen_events;
  float pot;
  float detector;
  float mc_type;
};

struct Slice {
  float flash_score;
  float flash_time;
};

struct FlatInteraction {
  PrimaryTrack ptrack;
  Slice slice;
  TrueNeutrino neutrino;
  EventInfo event;
  EventMeta meta;
};

} // end namespace
#endif
