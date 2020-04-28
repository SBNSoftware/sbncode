#include <iostream>

#include "PMTTrigger.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.h"

std::vector<numu::FlashTriggerPrimitive> numu::TriggerPrimitives(const std::vector<raw::OpDetWaveform> &waveforms, double tick_period, std::pair<double, double> &window, int thresh, bool is_sbnd) {
  std::vector<numu::FlashTriggerPrimitive> ret;

  opdet::sbndPDMapAlg mapalg;
  for (const raw::OpDetWaveform &wvf: waveforms) {
    // check if this is a PMT
    bool is_pmt = (!is_sbnd) || mapalg.isPDType(wvf.ChannelNumber(), "pmt");

    if (!is_pmt) continue;

    // look for any overlap with the enable window
    double waveform_start = wvf.TimeStamp();
    int waveform_index_start = std::max((int)((window.first - waveform_start) / tick_period), 0);
    int waveform_index_end = std::min((int)((window.second - waveform_start) / tick_period), (int) wvf.size());

    if (waveform_index_start < waveform_index_end) {
      numu::FlashTriggerPrimitive prim;
      prim.channel = wvf.ChannelNumber();
      for (int i = waveform_index_start; i < waveform_index_end; i++) {
        if (wvf[i] <= thresh) { // PMT waveforms go down
          FlashTriggerPrimitive::Trig this_trig {wvf[i], i + (int)((waveform_start - window.first) / tick_period)};
          prim.triggers.push_back(this_trig);
        }
      }
      ret.push_back(prim);
    }
  }
  return ret;
}

std::vector<int> numu::TriggerThresholds(const std::vector<numu::FlashTriggerPrimitive> &primitives, unsigned size) {
  std::map<int, std::vector<int>> trigs;
  for (const numu::FlashTriggerPrimitive &primitive: primitives) {
    for (const numu::FlashTriggerPrimitive::Trig &trig: primitive.triggers) {
      trigs[trig.tdc].push_back(trig.adc);
    }
  }

  std::map<int, std::pair<unsigned, int>> time_sorted;
  for (const auto &pair: trigs) {
    time_sorted[pair.first].first = pair.second.size();
    time_sorted[pair.first].second = *std::max_element(pair.second.begin(), pair.second.end());
  }

  std::vector<int> ret;
  for (unsigned i = 1; i <= size; i++) {
    int min = 8000;
    for (const auto &pair: time_sorted) {
      if (pair.second.first >= i) {
        if (pair.second.second < min) min = pair.second.second;
      }
    } 
    ret.push_back(min); 
  }

  return ret;
}

bool numu::HasTrigger(const std::vector<numu::FlashTriggerPrimitive> &primitives, int threshold, unsigned n_above_threshold) {
  if (n_above_threshold == 0) return true;

  std::map<int, std::vector<unsigned>> above_threshold;

  for (const numu::FlashTriggerPrimitive &primitive: primitives) {
    for (const numu::FlashTriggerPrimitive::Trig &trig: primitive.triggers) {
      if (trig.adc <= threshold) {
        above_threshold[trig.tdc].push_back(primitive.channel);
      }
    }
  }

  for (auto const &pair: above_threshold) {
    std::cout << "At time: " << pair.first << " above thresh: ";
    for (unsigned ch: pair.second) std::cout << ch << " ";
    std::cout << std::endl;
    if (pair.second.size() >= n_above_threshold) {
      return true;
    }
  }

  return false;
}


