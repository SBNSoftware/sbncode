/**
 * \file CRTSimHitViewer.cc
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

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTChannelData.h"
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "fhiclcpp/ParameterSet.h"

#include "TH2D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class CRTSimHitViewer
 * \brief Electron neutrino event selection
 */
class CRTSimHitViewer : public core::SelectionBase {
public:
  /** Constructor. */
  CRTSimHitViewer() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fCRTHitTag = config ? config->get<std::string>("CRTHitTag", "crtsimhit") : "crtsimhit";
    fIsICARUS = config ? config->get<bool>("IsICARUS", true) : true;
    event_ind = 0;
    
    time_v_pe_all = new TH2D("in_time_hits_all", "in_time_hits_all", 200, -1000., 1000., 10, 0., 500000.);
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
    fOutputFile->cd();
    time_v_pe_all->Write();
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {

    if (fIsICARUS) {
    auto const &crt_hits_handle = ev.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitTag);;
    const std::vector<sbn::crt::CRTHit> &crt_hits = *crt_hits_handle;
    const std::vector<sim::AuxDetSimChannel> &crt_sim_channels = *ev.getValidHandle<std::vector<sim::AuxDetSimChannel>>("largeant");
    art::FindManyP<icarus::crt::CRTData, void> hits_to_data(crt_hits_handle, ev, fCRTHitTag);
    std::cout << "Event: " << event_ind << std::endl;

    for (const sim::AuxDetSimChannel &sim_channel: crt_sim_channels) {
      int channel = sim_channel.AuxDetID();
      std::cout << "Channel: " << channel << std::endl;
      for (const sim::AuxDetIDE &ide: sim_channel.AuxDetIDEs()) {
        std::cout << "IDE\n";
        std::cout << "entry T: " << ide.entryT << " exit T: " << ide.exitT << std::endl; 
        std::cout << "entry X: " << ide.entryX << " entry Y: " << ide.entryY << " entry Z: " << ide.entryZ << std::endl;
        std::cout << "exit X: " << ide.exitX << " exit Y: " << ide.exitY << " exit Z: " << ide.exitZ << std::endl;
        std::cout << "Track ID: " << ide.trackID << std::endl;
      }
    }

    std::string histo_name = "in_time_hits" + std::to_string(event_ind);
    fOutputFile->cd();
    TH2D *time_v_pe = new TH2D(histo_name.c_str(), "in_time_hits", 200, -1000., 1000., 10, 0., 500000.);
    for (unsigned i = 0; i < crt_hits.size(); i++) {
      const sbn::crt::CRTHit &crt_hit = crt_hits[i];
      std::cout << "Hit index: " << i << std::endl;
      std::cout << "TS0 NS: " << (int)crt_hit.ts0_ns << std::endl;
      std::cout << "TS1 NS: " << crt_hit.ts1_ns << std::endl;
      std::cout << "X: " << crt_hit.x_pos << " Y: " << crt_hit.y_pos << " Z: " << crt_hit.z_pos << std::endl;
      std::cout << "Xerr: " << crt_hit.x_err << " Yerr: " << crt_hit.y_err << " Z: " << crt_hit.z_err << std::endl;
      std::cout << "Plane: " << crt_hit.plane << std::endl;
      std::cout << "PE's Hit: " << crt_hit.peshit << std::endl;
      std::cout << "Tagger: " << crt_hit.tagger << std::endl;
      time_v_pe->Fill(((int)crt_hit.ts0_ns) / 1000. + 1.1e3, crt_hit.peshit);
      time_v_pe_all->Fill(((int)crt_hit.ts0_ns) / 1000. + 1.1e3, crt_hit.peshit);
      const std::vector<art::Ptr<icarus::crt::CRTData>> &crt_datas = hits_to_data.at(i);
      std::cout << "Track IDs: ";
      for (auto const &crt_data: crt_datas) {
        std::vector<icarus::crt::CRTChannelData> crt_channel_datas = crt_data->ChanData();
        for (auto const &crt_chan: crt_channel_datas) {
          for (int id: crt_chan.TrackID()) {
            std::cout << id << " ";
          }
        }
      }
      std::cout << std::endl;
    }
    time_v_pe->Write();
    delete time_v_pe;
    }
    else {
    auto const &crt_hits_handle = ev.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitTag);;
    const std::vector<sbn::crt::CRTHit> &crt_hits = *crt_hits_handle;
    const std::vector<sim::AuxDetSimChannel> &crt_sim_channels = *ev.getValidHandle<std::vector<sim::AuxDetSimChannel>>("largeant");
    std::cout << "Event: " << event_ind << std::endl;

    for (const sim::AuxDetSimChannel &sim_channel: crt_sim_channels) {
      int channel = sim_channel.AuxDetID();
      std::cout << "Channel: " << channel << std::endl;
      for (const sim::AuxDetIDE &ide: sim_channel.AuxDetIDEs()) {
        std::cout << "IDE\n";
        std::cout << "entry T: " << ide.entryT << " exit T: " << ide.exitT << std::endl; 
        std::cout << "entry X: " << ide.entryX << " entry Y: " << ide.entryY << " entry Z: " << ide.entryZ << std::endl;
        std::cout << "exit X: " << ide.exitX << " exit Y: " << ide.exitY << " exit Z: " << ide.exitZ << std::endl;
      }
    }

    for (unsigned i = 0; i < crt_hits.size(); i++) {
      const sbn::crt::CRTHit &crt_hit = crt_hits[i];
      std::cout << "Hit index: " << i << std::endl;
      std::cout << "TS0 NS: " << crt_hit.ts0_ns << std::endl;
      std::cout << "TS1 NS: " << crt_hit.ts1_ns << std::endl;
      std::cout << "X: " << crt_hit.x_pos << " Y: " << crt_hit.y_pos << " Z: " << crt_hit.z_pos << std::endl;
      std::cout << "Xerr: " << crt_hit.x_err << " Yerr: " << crt_hit.y_err << " Z: " << crt_hit.z_err << std::endl;
      std::cout << "Plane: " << crt_hit.plane << std::endl;
      std::cout << "PE's Hit: " << crt_hit.peshit << std::endl;
    }

    }
    event_ind += 1;

    return false; 
  }

protected:
  std::string fCRTHitTag;
  bool fIsICARUS;
  unsigned event_ind;
  TH2D *time_v_pe_all;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::CRTSimHitViewer)

