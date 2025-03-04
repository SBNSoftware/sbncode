#include <iostream>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/ProductRetriever.h"
#include "art/Persistency/Common/GroupQueryResult.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"

#include "SBNDPOTTools.h"

using namespace std;

namespace sbn{
  sbn::PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags, int HLT) {
    bool foundHLT = false;
    PTBInfo_t PTBInfo;
    for (auto const& cont : *cont_frags)
    {
      artdaq::ContainerFragment cont_frag(cont);
      for (size_t fragi = 0; fragi < cont_frag.block_count(); ++fragi)
      {
        artdaq::Fragment frag = *cont_frag[fragi];
        sbndaq::CTBFragment ctb_frag(frag);   // somehow the name CTBFragment stuck
        for(size_t word_i = 0; word_i < ctb_frag.NWords(); ++word_i)
        {
          if(ctb_frag.Trigger(word_i)){
            if (ctb_frag.Trigger(word_i)->IsHLT() && ctb_frag.Trigger(word_i)->IsTrigger(HLT))
            {
              foundHLT = true;
              uint64_t RawprevPTBTimeStamp = ctb_frag.PTBWord(word_i)->prevTS * 20;
              uint64_t RawcurrPTBTimeStamp = ctb_frag.Trigger(word_i)->timestamp * 20;
              double currTS_candidate = std::bitset<64>(RawcurrPTBTimeStamp/20).to_ullong()/50e6;
              if(currTS_candidate < PTBInfo.currPTBTimeStamp){
                PTBInfo.prevPTBTimeStamp = std::bitset<64>(RawprevPTBTimeStamp / 20).to_ullong()/50e6;
                PTBInfo.currPTBTimeStamp = currTS_candidate;
                PTBInfo.GateCounter = ctb_frag.Trigger(word_i)->gate_counter;
              }
            }
          }
        } //End of loop over the number of trigger words
      } //End of loop over the number of fragments per container
    } //End of loop over the number of containers
  
    if(foundHLT == true){
      return PTBInfo;
    }
    else{
      std::cout << "Failed to find HLT 2!" << std::endl;
      throw std::exception();
    }
  }

  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) {
  
    double TDCTimeStamp = 0;
    for (auto const& cont : *cont_frags)
    {
      artdaq::ContainerFragment cont_frag(cont);
      for (size_t fragi = 0; fragi < cont_frag.block_count(); ++fragi)
      {
        artdaq::Fragment frag = *cont_frag[fragi];
        sbndaq::TDCTimestampFragment tdc_frag(frag);
        TDCTimeStamp = static_cast<double>(tdc_frag.getTDCTimestamp()->timestamp_ns())/1e9;
      } //End of loop over the number of fragments per container
    } //End of loop over the number of containers
    return TDCTimeStamp;
  }

}
