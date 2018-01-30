#include <cassert>
#include <json/json.h>
#include "Config.h"

namespace ana {
  namespace lee_truth_selection {

Config::Config(Json::Value* config) {
  Initialize(config);
}


void Config::Initialize(Json::Value* config) {
  Json::Value cfg;

  if (config) {
    cfg = *config;
  }

  // Global config
  ntrials = cfg.get("ntrials", 1).asInt();
  dataset_id = cfg.get("dataset_id", 0).asInt();
  track_energy_distortion = cfg.get("track_energy_distortion", 0).asFloat();
  shower_energy_distortion = cfg.get("shower_energy_distortion", 0).asFloat();
  do_misid = cfg.get("do_misid", false).asBool();

  // Selections
  accept_1l1p = cfg.get("accept_1l1p", true).asBool();
  accept_1l0p = cfg.get("accept_1l0p", false).asBool();
  accept_1lnp = cfg.get("accept_1lnp", false).asBool();
  accept_1lntrk = cfg.get("accept_1lntrk", false).asBool();

  if (accept_1l0p) {
    selections.push_back(k0p);
  }
  if (accept_1l1p) {
    selections.push_back(k1p);
  }
  if (accept_1lnp) {
    selections.push_back(kNp);
  }
  if (accept_1lntrk) {
    selections.push_back(kNtrk);
  }

  // Producers
  event_weight_producer = cfg.get("event_weight_producer", "").asString();
  mctruth_producer = cfg.get("mctruth_producer", "generator").asString();
  mcshower_producer = cfg.get("mcshower_producer", "mcreco").asString();
  mctrack_producer = cfg.get("mctrack_producer", "mcreco").asString();

  /*
  bool has_misids = false;
  std::vector<float> energy_range;
  std::vector<std::vector<int> > true_pdgids;
  std::vector<std::vector<int> > test_pdgids;
  std::vector<std::vector<float> > id_rates;
  if (cfg->isMember("do_misid")) {
    has_misids = true;
    auto inner = (*cfg)["do_misid"];
    for (auto energy: inner["energy_range"]) {
      float value = energy.asFloat();
      energy_range.push_back(value);

      true_pdgids.push_back(std::vector<int>());
      test_pdgids.push_back(std::vector<int>());
      id_rates.push_back(std::vector<float>());
      int size = id_rates.size();
      for (auto id_tuple: inner["id_rates"][energy.asCString()]) {
        true_pdgids[size-1].push_back(id_tuple[0].asInt());
        test_pdgids[size-1].push_back(id_tuple[1].asInt());
        id_rates[size-1].push_back(id_tuple[2].asFloat());
      }
    } 
  }*/
}

  }  // namespace lee_truth_selection
}  // namespace ana

