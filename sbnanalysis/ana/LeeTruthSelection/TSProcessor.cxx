#include <json/json.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TKey.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TTree.h>

#include "core/SelectionBase.hh"

#include "TSUtil.h"
#include "TSSelection.h"
#include "TSProcessor.h"

using namespace std;
using namespace art;

namespace ana {
namespace lee_truth_selection {

void TSProcessor::Initialize(Json::Value *config) {
  // All of the options
  int n_selections = config->get("N Selections", 1).asUInt();
  int n_trials = config->get("N Trials", 1).asUInt();
  
  std::vector<int> dataset_id;
  if (config->isMember("Dataset Id")) {
    for (auto id: (*config)["Datast Id"]) {
      dataset_id.push_back(id.asInt());
    }
  }

  std::vector<float> track_energy_distortion;
  if (config->isMember("Track Energy Distortion")) {
    for (auto dist: (*config)["Track Energy Distortion"]) {
      track_energy_distortion.push_back(dist.asFloat());
    }
  }
  bool track_energy_distortion_by_percent = 
    config->get("Track Energy Distortion By Percent", true).asBool();

  std::vector<float> shower_energy_distortion;
  if (config->isMember("Shower Energy Distortion")) {
    for (auto dist: (*config)["Shower Energy Distortion"]) {
      shower_energy_distortion.push_back(dist.asFloat());
    }
  }
  bool shower_energy_distortion_by_percent = 
    config->get("Shower Energy Distortion By Percent", true).asBool();

  bool drop_np = config->get("Drop NP", false).asBool();
  bool drop_ntrack = config->get("Drop NTrk", false).asBool();

  bool has_misids = false;
  std::vector<float> energy_range;
  std::vector<std::vector<int>> true_pdgids;
  std::vector<std::vector<int>> test_pdgids;
  std::vector<std::vector<float>> id_rates;
  if (config->isMember("Mis-id")) {
    has_misids = true;
    auto inner = (*config)["Mis-id"];
    for (auto energy: inner["Energy Range"]) {
      float value = energy.asFloat();
      energy_range.push_back(value);

      true_pdgids.push_back(std::vector<int>());
      test_pdgids.push_back(std::vector<int>());
      id_rates.push_back(std::vector<float>());
      int size = id_rates.size();
      for (auto id_tuple: inner["ID Rates"][energy.asCString()]) {
        true_pdgids[size-1].push_back(id_tuple[0].asInt());
        test_pdgids[size-1].push_back(id_tuple[1].asInt());
        id_rates[size-1].push_back(id_tuple[2].asFloat());
      }
    } 
  }

  std::vector<TFile *> f_outs;
  if (config->isMember("Output Files")) {
    for (auto f_name: (*config)["Output Files"]) {
      f_outs.push_back(new TFile(f_name.asCString()));
    }
  }

  cout << "Initialize" << endl;
  for (int i = 0; i < n_selections; i ++) {
    _selections[i].setOutputFile(f_outs[i]); 
    // TODO: populate this vector
    std::vector<std::string> filenames;
    _selections[i].initialize(filenames);
  }

  for (int i = 0; i < n_selections; i ++) {
    _selections[i].setFluxWeightProducer("eventweight");
    _selections[i].setEventWeightProducer("mcweight");
    _selections[i].setMCTruthProducer("generator");
    _selections[i].setMCTrackProducer("mcreco");
    _selections[i].setMCShowerProducer("mcreco");
    _selections[i].setVerbose(true);
    _selections[i].setNTrials(n_trials);
    _selections[i].setAcceptP(!drop_np, 2);
    _selections[i].setAcceptNTrk(!drop_ntrack);

    if (has_misids) {
      _selections[i].setParticleIDEnergyRange(energy_range);
      for (unsigned j = 0; j < energy_range.size(); j++) {
        for (unsigned k = 0; k < id_rates[j].size(); k++) {
          int true_pdgid = true_pdgids[j][k];
          int test_pdgid = test_pdgids[j][k];
          float id_rate = id_rates[j][k];
          _selections[i].addParticleIDRate(true_pdgid, test_pdgid, id_rate, energy_range[j] - 0.1);
          cout << "At ENERGY: " << energy_range[j] << " TRUE ID " << true_pdgid << " TEST ID " << test_pdgid << " RATE " << id_rate << endl;
        }
      }
      _selections[i].checkParticleIDRates();
    }
   
    if (dataset_id.size() > 0) {
        _selections[i].setDatasetID( dataset_id[i] );
    } 
    if (track_energy_distortion.size() > 0) {
      _selections[i].setTrackEnergyResolution(track_energy_distortion[i], track_energy_distortion_by_percent);
    } 
    if (shower_energy_distortion.size() > 0) {
      _selections[i].setShowerEnergyResolution(shower_energy_distortion[i], shower_energy_distortion_by_percent);
    }
  }
}

void TSProcessor::ProcessEvent(gallery::Event& ev) {
  for (ana::lee_truth_selection::TSSelection &selection: _selections) {
    selection.analyze(&ev);
  }
}
 
void TSProcessor::Finalize() { 
  for (ana::lee_truth_selection::TSSelection &selection: _selections) {
    selection.finalize();
  }
}

}
}
