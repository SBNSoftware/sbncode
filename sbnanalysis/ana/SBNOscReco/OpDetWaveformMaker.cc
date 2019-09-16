/**
 * \file OpDetWaveformMaker.cc
 *
 *
 * Author:
 */

#include <iostream>
#include <array>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/Experiment.hh"
#include "core/ProviderManager.hh"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class OpDetWaveformMaker
 * \brief Electron neutrino event selection
 */
class OpDetWaveformMaker : public core::SelectionBase {
public:
  /** Constructor. */
  OpDetWaveformMaker() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fOpDetWaveformTag = config ? config->get<std::string>("OpDetWaveformTag", "opdaq") : "opdaq";
    clock_freq = fProviderManager->GetDetectorClocksProvider()->OpticalClock().Frequency() * (500. / 64 /* bug fix???*/); // / 1000. /* convert to GHz */;
    std::cout << "CLOCK FREQ: " << clock_freq << std::endl;
    event_ind = 0;
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
    fOutputFile->cd();
    for (TGraph *plot: plots) {
      plot->Write();
    }
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
    const std::vector<raw::OpDetWaveform> &waveforms = *ev.getValidHandle<std::vector<raw::OpDetWaveform>>(fOpDetWaveformTag);
    //std::map<raw::Channel_t, std::pair<std::vector<double>, std::vector<double>>> waveform_plots;
    std::map<raw::Channel_t, std::vector<std::pair<double, double>>> waveform_data;
    std::map<raw::Channel_t, std::vector<double>> wf_start_times;
    std::map<raw::Channel_t, std::vector<double>> wf_periods;
    for (const raw::OpDetWaveform &wfconst: waveforms) {
      raw::OpDetWaveform wf = wfconst;
      //if (waveform_plots.count(wf.ChannelNumber()) > 0) {
      //  std::cout << "Multi-waveform output! Channel: " << wf.ChannelNumber() << " Event: " << event_ind << std::endl;
      //}
      for (unsigned i = 0; i < wf.Waveform().size(); i++) {
        waveform_data[wf.ChannelNumber()].push_back({i / clock_freq + wf.TimeStamp(), (double)wf.Waveform()[i]});
      }
      wf_start_times[wf.ChannelNumber()].push_back(wf.TimeStamp());
      wf_periods[wf.ChannelNumber()].push_back( wf.Waveform().size() / clock_freq);
    }

    // now sort each one
    //for (auto &wf_data_pair: waveform_data) {
    //  std::sort(wf_data_pair.second.begin(), wf_data_pair.second.end(), [](auto lhs, auto rhs) {return lhs.first < rhs.first;});
    //}
    std::map<raw::Channel_t, std::pair<std::vector<double>, std::vector<double>>> waveform_plots;
    for (auto &wf_data_pair: waveform_data) {
      for (auto data_pair: wf_data_pair.second) {
        waveform_plots[wf_data_pair.first].first.push_back(data_pair.first);
        waveform_plots[wf_data_pair.first].second.push_back(data_pair.second);
      }
    }

    for (auto &start_time_pair: wf_start_times) {
      std::cout << "Event: " << event_ind << " Channel: " << start_time_pair.first << std::endl;
      const std::vector<double> &wf_period = wf_periods.at(start_time_pair.first); 
      for (unsigned i = 0; i < start_time_pair.second.size(); i++) {
        std::cout << start_time_pair.second[i] << " (" << wf_period[i] <<") "; 

      }
      std::cout << std::endl;
    }

    for (auto &plot_info: waveform_plots) {
      TGraph *plot = new TGraph(plot_info.second.first.size(), &plot_info.second.first[0], &plot_info.second.second[0]);
      plot->SetName((std::string("Waveform ch: ") + std::to_string(plot_info.first) + " ev: " + std::to_string(event_ind)).c_str());
      plot->SetTitle((std::string("Waveform ch: ") + std::to_string(plot_info.first) + " ev: " + std::to_string(event_ind)).c_str());
      plots.push_back(plot);
    }
    event_ind += 1;

    return false; 
  }

protected:
  std::string fOpDetWaveformTag;
  double clock_freq;
  unsigned event_ind;
  std::vector<TGraph *> plots;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::OpDetWaveformMaker)

