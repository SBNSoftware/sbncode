/**
 * \file CRTHitTiming.cc
 *
 *
 * Author:
 */

#include <iostream>
#include <array>

#include "nusimdata/SimulationBase/MCParticle.h"

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

#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"

#include "TH2D.h"
#include "TGraph.h"

namespace ana {
  namespace SBNOsc {

/**
 * \class CRTHitTiming
 * \brief Electron neutrino event selection
 */
class CRTHitTiming : public core::SelectionBase {
public:
  /** Constructor. */
  CRTHitTiming() {}

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL) {
    fCRTHitTag = config ? config->get<std::string>("CRTHitTag", "crtsimhit") : "crtsimhit";
    fIsICARUS = config ? config->get<bool>("IsICARUS", true) : true;
    event_ind = 0;
    _verbose = false;
    _true_v_hit_time = new TH1D("true_v_hit_time", "true_v_hit_time", 200, -100., 100.);
    _true_v_hit_time_c = new TH1D("true_v_hit_time_c", "true_v_hit_time_c", 200, -100., 100.);
    _true_v_hit_time_d = new TH1D("true_v_hit_time_d", "true_v_hit_time_d", 200, -100., 100.);
    _true_v_hit_time_m = new TH1D("true_v_hit_time_m", "true_v_hit_time_m", 200, -100., 100.);
    _true_v_hit_time_coince = new TH1D("true_v_hit_time_coince", "true_v_hit_time_coince", 200, -100., 100.);
    _hitpe_muon = new TH1D("hitpe_muon", "hitpe_muon", 1000, 0., 1000.);
    _hitpe_proton = new TH1D("hitpe_proton", "hitpe_proton", 1000, 0., 1000.);
    _hitpe_neutron = new TH1D("hitpe_neutron", "hitpe_neutron", 1000, 0., 1000.);
    FillFebMap();
  }

  /** Finalize and write objects to the output file. */
  void Finalize() {
    fOutputFile->cd();
    _true_v_hit_time_coince->Write();
    _true_v_hit_time->Write();
    _true_v_hit_time_c->Write();
    _true_v_hit_time_m->Write();
    _true_v_hit_time_d->Write();
    _hitpe_muon->Write();
    _hitpe_proton->Write();
    _hitpe_neutron->Write();
  }

  int GetPDG(const std::vector<simb::MCParticle> &mc_particles, int trackID) {
    for (const simb::MCParticle &part: mc_particles) {
      if (part.TrackId() == trackID) return part.PdgCode();
    }
    return -1;
  }

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
      const std::vector<simb::MCParticle> &mc_particles = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");

    if (fIsICARUS) {
      auto const &crt_hits_handle = ev.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitTag);;
      const std::vector<sbn::crt::CRTHit> &crt_hits = *crt_hits_handle;
      const std::vector<sim::AuxDetSimChannel> &crt_sim_channels = *ev.getValidHandle<std::vector<sim::AuxDetSimChannel>>("largeant");
      art::FindManyP<icarus::crt::CRTData, void> hits_to_data(crt_hits_handle, ev, fCRTHitTag);

      // TODO: fix when ICARUS has the necessary associations
      // collect the mapping of track ID to IDE's
      std::map<int, std::vector<std::pair<std::pair<int, int>, sim::AuxDetIDE>>> tracks_to_ides;
      for (const sim::AuxDetSimChannel &simchannel: crt_sim_channels) {
        for (const sim::AuxDetIDE &ide: simchannel.AuxDetIDEs()) {
          tracks_to_ides[ide.trackID].push_back({{simchannel.AuxDetID(), simchannel.AuxDetSensitiveID()}, ide});
        }
      }

      // Now for each CRT hit, try to figure out the IDE(s) that made it
      for (unsigned i = 0; i < crt_hits.size(); i++) {
        const sbn::crt::CRTHit &hit = crt_hits[i];
        std::vector<art::Ptr<icarus::crt::CRTData>> datas = hits_to_data.at(i);
        assert(datas.size() == 1 || datas.size() == 2);
        bool has_coincedence = datas.size() == 2;
        _verbose = false;

        // map to auxdet
        int auxdetid_0 = MacToAuxDetID(datas[0]->Mac5(), datas[0]->ChanTrig());
        // and to type
        const geo::AuxDetGeo& adGeo = fProviderManager->GetAuxDetGeometryProvider()->AuxDet(auxdetid_0);
        const char auxDetType = GetAuxDetType(adGeo); //CRT module type (c, d, or m)

        if (has_coincedence) {
          int auxdetid_1 = MacToAuxDetID(datas[1]->Mac5(), datas[1]->ChanTrig());
          const geo::AuxDetGeo& adGeo = fProviderManager->GetAuxDetGeometryProvider()->AuxDet(auxdetid_1);
          assert(auxDetType == GetAuxDetType(adGeo));
        }

        double hit_time = (int)hit.ts0_ns;

        if (_verbose) std::cout << "CRT Hit at: " << hit.ts1_ns << " " << (int)hit.ts0_ns << " " << (int)hit.ts0_ns_corr << std::endl;
        if (_verbose) std::cout << "Has Coincedence: " << has_coincedence << std::endl;

        double true_time_sum = 0.;
        unsigned n_true_time = 0;

        bool has_muon = false;
        bool has_proton = false;
        bool has_neutron = false;
        for (unsigned j = 0; j < datas.size(); j++) {
          const icarus::crt::CRTData &data = *datas[j];
          if (_verbose) std::cout << "CRT Data at: " << data.TTrig() << " Mac5: " << data.Mac5() << " chan: " << data.ChanTrig() << std::endl;
          const std::vector<icarus::crt::CRTChannelData> &channel_datas = data.ChanData();

        
          for (const icarus::crt::CRTChannelData &chan: channel_datas) {
            const std::vector<int> track_ids = chan.TrackID();
            if (_verbose) std::cout << "Channel Data at: " << chan.T0() << " chan: " << chan.Channel() << std::endl; 
            for (int track_id: track_ids) {
              int pdg_match = abs(GetPDG(mc_particles, track_id));
              if (abs(pdg_match) == 13) has_muon = true;
              else if (abs(pdg_match) == 2212) has_proton = true;
              else if (abs(pdg_match) == 2112) has_neutron = true;
              if (_verbose) std::cout << "True Particle: " << track_id << std::endl;
              for (const auto &ide_pair: tracks_to_ides.at(track_id)) {
                const sim::AuxDetIDE &ide = ide_pair.second;
                std::pair<int, int> channel = ide_pair.first;
                double true_time = ide.entryT - 1.1e6;
                if (_verbose) std::cout << "True Hit at: " << true_time << " X: "  << ide.entryX << " Y: " << ide.entryY <<" Z: " << ide.entryZ << " Chan: " << channel.first << " " << channel.second << std::endl;
                std::array<int, 3> crt_channel = GetCRTChannel(channel.first, channel.second);
                if (_verbose) std::cout << "Mapped FEB: " << crt_channel[0] << "Mapped channel: " << crt_channel[1] << " " << crt_channel[2] << std::endl;

                if (data.Mac5() == crt_channel[0] && (data.ChanTrig() == crt_channel[1] || data.ChanTrig() == crt_channel[2])) {
                  true_time_sum += true_time;
                  n_true_time += 1;
                }
              }
            }
          }
        }
        if (has_coincedence) _true_v_hit_time_coince->Fill((true_time_sum/n_true_time) - hit_time);
        else {
          std::cout << "Hit time: " << hit_time << " avg true time: " << (true_time_sum/n_true_time) << std::endl;
          _true_v_hit_time->Fill((true_time_sum/n_true_time) - hit_time);
        }
        if (auxDetType == 'm') _true_v_hit_time_m->Fill((true_time_sum/n_true_time) - hit_time);
        else if (auxDetType == 'd') _true_v_hit_time_d->Fill((true_time_sum/n_true_time) - hit_time);
        else if (auxDetType == 'c') _true_v_hit_time_c->Fill((true_time_sum/n_true_time) - hit_time);
        else assert(false);

        double pes = (hit.peshit - 2*63.6) / (70.*2);
        if (has_muon) _hitpe_muon->Fill(pes);
        else if (has_proton) _hitpe_proton->Fill(pes);
        else if (has_neutron) _hitpe_neutron->Fill(pes);
        
      }
      // fOutputFile->cd();
    }
    else {
      auto const &crt_hits_handle = ev.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitTag);;
      const std::vector<sbn::crt::CRTHit> &crt_hits = *crt_hits_handle;
      const std::vector<sim::AuxDetSimChannel> &crt_sim_channels = *ev.getValidHandle<std::vector<sim::AuxDetSimChannel>>("largeant");
      art::FindManyP<sbnd::crt::CRTData, void> hits_to_data(crt_hits_handle, ev, fCRTHitTag);
      auto const &crt_data_handle = ev.getValidHandle<std::vector<sbnd::crt::CRTData>>("crt");
      for (unsigned i = 0; i < crt_hits.size(); i++) {
        int track_id = -1;
        const sbn::crt::CRTHit &hit = crt_hits[i];
        std::vector<art::Ptr<sbnd::crt::CRTData>> datas = hits_to_data.at(i);
        art::FindManyP<sim::AuxDetIDE, void> sims(datas, ev, "crt");
        bool has_muon = false;
        bool has_proton = false;
        bool has_neutron = false;
        for (unsigned j = 0; j < datas.size(); j++) {
          for (unsigned k = 0; k < sims.at(j).size(); k++) {
            track_id = sims.at(j).at(k)->trackID;
            int pdg_match = GetPDG(mc_particles, track_id);
	    if (abs(pdg_match) == 13) has_muon = true;
	    else if (abs(pdg_match) == 2212) has_proton = true;
	    else if (abs(pdg_match) == 2112) has_neutron = true;
          }
        }
        double pes = hit.peshit;
        if (has_muon) _hitpe_muon->Fill(pes);
        else if (has_proton) _hitpe_proton->Fill(pes);
        else if (has_neutron) _hitpe_neutron->Fill(pes);
      }
    }
    event_ind += 1;

    return false; 
  }

protected:
  void FillFebMap()
  {
    if(!this->fFebMap.empty())
        return;

    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    searchPath.find_file("feb_map.txt",fullFileName);
    std::ifstream fin;
    fin.open(fullFileName,std::ios::in);
    if(fin.good()) std::cout << "opened file 'feb_map.txt' for reading..." << std::endl;
    else
        throw cet::exception("FillFebMap") << "Unable to find/open file 'feb_map.txt'" << std::endl;
    std::vector<std::string> row;
    std::string line, word;
    while(getline(fin,line)) {
        row.clear();
        std::stringstream s(line);
        int mod;
        while (std::getline(s, word, ',')) {
            row.push_back(word);
        }
        mod = std::stoi(row[0]);
        (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[1]),std::stoi(row[2])));
        if(row.size()>3)
            (this->fFebMap)[mod].push_back(std::make_pair(std::stoi(row[3]),std::stoi(row[4])));
    }
    std::cout << "filled febMap with " << (this->fFebMap).size() << " entries" << std::endl;
    fin.close();
  } //FillFebMap()

char GetAuxDetType(geo::AuxDetGeo const& adgeo)
{
  std::string volName(adgeo.TotalVolume()->GetName());
  if (volName.find("MINOS") != std::string::npos) return 'm';
  if (volName.find("CERN")  != std::string::npos) return 'c';
  if (volName.find("DC")    != std::string::npos) return 'd';

  return 'e';
}//GetAuxDetType()

  char MacToType(int mac)
  {

      int reg = MacToRegion(mac);
      if(reg>=30&&reg<40) return 'c';
      if(reg>=40&&reg<50) return 'm';
      if(reg==50) return 'd';
      std::cout << "ERROR in MacToType: type not set!" << std::endl;
      return 'e';
  }

  int MacToRegion(int mac){

      if(mac>=107 && mac<=190) return 30; //top
      if(mac>=191 && mac<=204) return 31; //rim west
      if(mac>=205 && mac<=218) return 32; //rim east
      if(mac>=219 && mac<=224) return 33; //rim south
      if(mac>=225 && mac<=230) return 34; //rim north
      if(            mac<=12 ) return 40; //west side, south stack
      if(mac>=13  && mac<=24 ) return 41; //west side, center stack
      if(mac>=25  && mac<=36 ) return 42; //west side, north stack
      if(mac>=37  && mac<=48 ) return 43; //east side, south stack
      if(mac>=49  && mac<=60 ) return 44; //east side, center stack
      if(mac>=61  && mac<=72 ) return 45; //east side, north stack
      if(mac>=73  && mac<=84 ) return 46; //south
      if(mac>=85  && mac<=92 ) return 47; //north
      if(mac>=93 && mac<=106) return 50; //bottom

      std::cout << "ERROR in MacToRegion: unknown mac address " << mac << std::endl;
      return 0;
  }

  int MacToAuxDetID(int mac, int chan)
  {
      char type = MacToType(mac);
      if (type == 'e') return INT_MAX;

      int pos=1;
      if(type=='m')
          pos = chan/10 + 1;

      if((this->fFebMap).empty())
          std::cout << "ERROR in MacToAuxDetID: FEBMap is empty!" << std::endl;

      for(const auto&  p : this->fFebMap) {
          if(p.second[0].first == mac && p.second[0].second==pos)
              return (uint32_t)p.first;
          if(p.second.size()==2)
              if(p.second[1].first==mac && p.second[1].second==pos)
                  return (uint32_t)p.first;
      }


    std::cout << "ERROR in MacToAuxDetID: auxDetID not set!" << std::endl;
    return INT_MAX;
  }


  std::array<int, 3> GetCRTChannel(int aux_det_id, int aux_det_sens_id) {
      // const char auxDetType, ) {
      const int adid = aux_det_id; //adsc.AuxDetID(); //CRT module ID number (from gdml)
      const int adsid = aux_det_sens_id; //adsc.AuxDetSensitiveID(); //CRT strip ID number (from gdml)
      const geo::AuxDetGeo& adGeo = fProviderManager->GetAuxDetGeometryProvider()->AuxDet(adid);
      const char auxDetType = GetAuxDetType(adGeo); //CRT module type (c, d, or m)

      int channel0ID=-1, channel1ID=-1;
      int mac5;
      int mac5dual;
      switch (auxDetType){
          case 'c' :
              mac5 =fFebMap[adid][0].first;
              channel0ID = 2 * adsid + 0;
              channel1ID = 2 * adsid + 1;
              if(mac5<107||mac5>230)
                  std::cout << "WARNING: mac5 out of bounds for c-type!" << std::endl;
              if(channel0ID<0 || channel0ID > 31 || channel1ID<0 || channel1ID>31)
                  std::cout << "WARNING: channel out of bounds for c-type!" << std::endl;
              break;
          case 'd' : 
              mac5 = fFebMap[adid][0].first;
              channel0ID = adsid;
              if(mac5<93||mac5>106)
                  std::cout << "WARNING: mac5 out of bounds for d-type!" << std::endl;
              if(channel0ID<0 || channel0ID > 63)
                  std::cout << "WARNING: channel out of bounds for d-type!" << std::endl;
              break;
          case 'm' : 
              mac5 = fFebMap[adid][0].first;
              channel0ID = adsid/2 + 10*(fFebMap[adid][0].second-1);
              if(mac5<1||mac5>92)
                  std::cout << "WARNING: mac5 out of bounds for m-type!" << std::endl;
              if(channel0ID<0 || channel0ID > 31)
                  std::cout << "WARNING: channel out of bounds for m-type!" << std::endl;
              if (fFebMap[adid].size()==2)  {
                  mac5dual = fFebMap[adid][1].first;
                  if(mac5dual<1||mac5dual>92)
                      std::cout << "WARNING: mac5dual out of bounds for m-type!" << std::endl;
              }
              break;
      }
      return {mac5, channel0ID, channel1ID};
  }

  std::map<int,std::vector<std::pair<int,int>>> fFebMap;
  std::string fCRTHitTag;
  bool fIsICARUS;
  bool _verbose;
  unsigned event_ind;
  TH1D *_true_v_hit_time;
  TH1D *_true_v_hit_time_coince;
  TH1D *_true_v_hit_time_c;
  TH1D *_true_v_hit_time_d;
  TH1D *_true_v_hit_time_m;
  TH1D *_hitpe_muon;
  TH1D *_hitpe_proton;
  TH1D *_hitpe_neutron;
};

  }  // namespace SBNOsc
}  // namespace ana
DECLARE_SBN_PROCESSOR(ana::SBNOsc::CRTHitTiming)

