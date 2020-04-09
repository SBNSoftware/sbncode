////////////////////////////////////////////////////////////////////////
// Class:       PMTFlashTriggerMaker
// Plugin Type: producer (art v3_02_06)
// File:        PMTFlashTriggerMaker_module.cc
//
// Generated at Wed Feb 19 17:38:21 2020 by Gray Putnam using cetskelgen
// from cetlib version v3_07_02.
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
#include "lardata/Utilities/AssociationUtil.h"

#include "sbncode/LArRecoProducer/Products/FlashTriggerPrimitive.hh"

#include <memory>
#include <iostream>
#include <vector>

#include "sbndcode/OpDetSim/sbndPDMapAlg.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace sbn {
  class PMTFlashTriggerMaker;

  std::vector<sbn::FlashTriggerPrimitive> TriggerPrimitives(
    const std::vector<raw::OpDetWaveform> &waveforms, 
    double tick_period, 
    std::pair<double, double> &window, 
    int thresh, 
    bool is_sbnd);
}


class sbn::PMTFlashTriggerMaker : public art::EDProducer {
public:
  explicit PMTFlashTriggerMaker(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTFlashTriggerMaker(PMTFlashTriggerMaker const&) = delete;
  PMTFlashTriggerMaker(PMTFlashTriggerMaker&&) = delete;
  PMTFlashTriggerMaker& operator=(PMTFlashTriggerMaker const&) = delete;
  PMTFlashTriggerMaker& operator=(PMTFlashTriggerMaker&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:
  art::InputTag fWaveformLabel;
  std::string fExperiment;
  std::pair<double, double> fTriggerWindow;
  int fTriggerThreshold;
  bool fOffsetTriggerTime;
};

sbn::PMTFlashTriggerMaker::PMTFlashTriggerMaker(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fWaveformLabel(p.get<std::string>("WaveformLabel")),
    fExperiment(p.get<std::string>("Experiment")),
    fTriggerWindow({p.get<float>("TriggerStart"), p.get<float>("TriggerEnd")}),
    fTriggerThreshold(p.get<int>("TriggerThreshold")),
    fOffsetTriggerTime(p.get<bool>("OffsetTriggerTime"))
{
  produces< std::vector<sbn::FlashTriggerPrimitive> >();
}

void sbn::PMTFlashTriggerMaker::produce(art::Event& e)
{
  std::unique_ptr<std::vector<sbn::FlashTriggerPrimitive>> trigs(new std::vector<sbn::FlashTriggerPrimitive>);

  art::Handle<std::vector<raw::OpDetWaveform>> waveform_handle;
  e.getByLabel(fWaveformLabel, waveform_handle);

  std::pair<double, double> window = fTriggerWindow;
  // don't apply offset in constructor in case time gets overwritten on input file open
  if (fOffsetTriggerTime) {
    float offset = lar::providerFrom<detinfo::DetectorClocksService>()->TriggerTime();
    window.first += offset;
    window.second += offset;
  }

  if (waveform_handle.isValid()) {
    const std::vector<raw::OpDetWaveform> &waveforms = *waveform_handle;
    double tick_period = lar::providerFrom<detinfo::DetectorClocksService>()->OpticalClock().TickPeriod();
    bool is_sbnd = fExperiment == "SBND";

   *trigs = TriggerPrimitives(waveforms, tick_period, window, fTriggerThreshold, is_sbnd);
  }

  e.put(std::move(trigs));
}

DEFINE_ART_MODULE(sbn::PMTFlashTriggerMaker)

std::vector<sbn::FlashTriggerPrimitive> sbn::TriggerPrimitives(const std::vector<raw::OpDetWaveform> &waveforms, double tick_period, std::pair<double, double> &window, int thresh, bool is_sbnd) {
  std::vector<sbn::FlashTriggerPrimitive> ret;

  opdet::sbndPDMapAlg mapalg;
  for (const raw::OpDetWaveform &wvf: waveforms) {
    // check if this is a PMT
    bool is_pmt = (!is_sbnd) || mapalg.pdType(wvf.ChannelNumber(), "pmt");

    if (!is_pmt) continue;

    // look for any overlap with the enable window
    double waveform_start = wvf.TimeStamp();
    int waveform_index_start = std::max((int)((window.first - waveform_start) / tick_period), 0);
    int waveform_index_end = std::min((int)((window.second - waveform_start) / tick_period), (int) wvf.size());

    if (waveform_index_start < waveform_index_end) {
      sbn::FlashTriggerPrimitive prim;
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
