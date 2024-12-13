////////////////////////////////////////////////////////////////////////
// Class:       SBNDBNBEXTRetriever
// Plugin Type: producer 
// File:        SBNDBNBEXTRetriever_module.cc
//
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

#include <memory>
#include <bitset>
#include <tuple>
#include <algorithm>

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
//#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"

#include "ifdh_art/IFBeamService/IFBeam_service.h"
#include "ifbeam_c.h"

#include "larcorealg/CoreUtils/counter.h"

namespace sbn {
  class SBNDBNBEXTRetriever;
}

class sbn::SBNDBNBEXTRetriever : public art::EDProducer {
public:
  explicit SBNDBNBEXTRetriever(fhicl::ParameterSet const & params);
  // Required functions.
  void produce(art::Event & e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

  // Plugins should not be copied or assigned.
  SBNDBNBEXTRetriever(SBNDBNBEXTRetriever const &) = delete;
  SBNDBNBEXTRetriever(SBNDBNBEXTRetriever &&) = delete;
  SBNDBNBEXTRetriever & operator = (SBNDBNBEXTRetriever const &) = delete;
  SBNDBNBEXTRetriever & operator = (SBNDBNBEXTRetriever &&) = delete;


private:
  // Declare member data here.
  std::vector< sbn::EXTCountInfo > fOutExtInfos;
  struct PTBInfo_t {
    double currPTBTimeStamp  = 0;
    double prevPTBTimeStamp  = 0;
    unsigned int GateCounter = 0; // FIXME needs to be integral type
  };

  struct TriggerInfo_t {
    double t_current_event  = 0;
    double t_previous_event = 0;
    unsigned int number_of_gates_since_previous_event = 0; // FIXME needs to be integral type
  };

  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  PTBInfo_t extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const;
  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const;

  // input labels
  std::string raw_data_label_;
  float TotalEXTCounts;  
  float totalMinBias;
  float evtCount;
};

sbn::SBNDBNBEXTRetriever::SBNDBNBEXTRetriever(fhicl::ParameterSet const & params)
  : EDProducer{params} {
  produces< std::vector< sbn::EXTCountInfo >, art::InSubRun >();
  TotalEXTCounts = 0;
  totalMinBias = 0;
  evtCount = 0;
}

int eventNum =0;
int _run;
int _subrun;
int _event;

void sbn::SBNDBNBEXTRetriever::produce(art::Event & e)
{


  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);

  TotalEXTCounts += triggerInfo.number_of_gates_since_previous_event;

  if(triggerInfo.number_of_gates_since_previous_event > 0){
    evtCount++;  
    totalMinBias += triggerInfo.number_of_gates_since_previous_event;
  }
   
  //Store everything in our data-product
  sbn::EXTCountInfo extInfo;
  extInfo.gates_since_last_trigger = triggerInfo.number_of_gates_since_previous_event;
  fOutExtInfos.push_back(extInfo);
}

sbn::SBNDBNBEXTRetriever::PTBInfo_t sbn::SBNDBNBEXTRetriever::extractPTBInfo(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const {
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
          uint32_t wt = 0;
          uint32_t word_type = ctb_frag.Word(word_i)->word_type;
          wt = word_type;
	  if (wt == 2 && ctb_frag.Trigger(word_i)->IsTrigger(4))
	  {
	    uint64_t RawprevPTBTimeStamp = ctb_frag.PTBWord(word_i)->prevTS * 20; 
	    uint64_t RawcurrPTBTimeStamp = ctb_frag.Trigger(word_i)->timestamp * 20;
            PTBInfo.prevPTBTimeStamp = std::bitset<64>(RawprevPTBTimeStamp / 20).to_ullong()/50e6; 
            PTBInfo.currPTBTimeStamp = std::bitset<64>(RawcurrPTBTimeStamp / 20).to_ullong()/50e6; 
            PTBInfo.GateCounter = ctb_frag.Trigger(word_i)->gate_counter;
            break;
	  }
        }
      } //End of loop over the number of trigger words
    } //End of loop over the number of fragments per container
  } //End of loop over the number of containers

  return PTBInfo; 
}

double sbn::SBNDBNBEXTRetriever::extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const {
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

sbn::SBNDBNBEXTRetriever::TriggerInfo_t sbn::SBNDBNBEXTRetriever::extractTriggerInfo(art::Event const& e) const {
  // Using TDC for current event, but PTB for previous event
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);

  art::InputTag TDC_itag("daq", "ContainerTDCTIMESTAMP");
  auto TDC_cont_frags = e.getHandle<artdaq::Fragments>(TDC_itag);

  PTBInfo_t PTBInfo;
  TriggerInfo_t triggerInfo;
  PTBInfo = extractPTBInfo(PTB_cont_frags);

  if (TDC_cont_frags) {
    double TDCTimeStamp = extractTDCTimeStamp(TDC_cont_frags);
    triggerInfo.t_current_event = TDCTimeStamp;
  }
  else{
    mf::LogDebug("SBNDBNBEXTRetriever") << " Missing TDC Contaienr Fragments!!! " << std::endl;
    triggerInfo.t_current_event = PTBInfo.currPTBTimeStamp;
  }

  triggerInfo.t_previous_event = PTBInfo.prevPTBTimeStamp;
  triggerInfo.number_of_gates_since_previous_event = PTBInfo.GateCounter;

  if(triggerInfo.t_current_event - PTBInfo.currPTBTimeStamp >= 1){
    mf::LogDebug("SBNDBNBEXTRetriever") << "Caught PTB bug, PTB late" << std::endl;
    mf::LogDebug("SBNDBNBEXTRetriever") << "Before: " << triggerInfo.t_previous_event << std::endl;
    triggerInfo.t_previous_event+=1;
    mf::LogDebug("SBNDBNBEXTRetriever") << "After: " << triggerInfo.t_previous_event << std::endl;
  }
  else if(triggerInfo.t_current_event - PTBInfo.currPTBTimeStamp <= -1){
    mf::LogDebug("SBNDBNBEXTRetriever") << "Caught PTB bug, PTB early" << std::endl;
    mf::LogDebug("SBNDBNBEXTRetriever") << "Before: " << triggerInfo.t_previous_event << std::endl;
    triggerInfo.t_previous_event-=1;
    mf::LogDebug("SBNDBNBEXTRetriever") << "After: " << triggerInfo.t_previous_event << std::endl;
  }

  return triggerInfo;
}

void sbn::SBNDBNBEXTRetriever::beginSubRun(art::SubRun& sr)
{
  TotalEXTCounts = 0;
  totalMinBias = 0;
  evtCount = 0;
  return;
}

void sbn::SBNDBNBEXTRetriever::endSubRun(art::SubRun& sr)
{
   // We will add all of the EXTCountInfo data-products to the 
  // art::SubRun so it persists 

  auto p =  std::make_unique< std::vector< sbn::EXTCountInfo > >(fOutExtInfos);

  sr.put(std::move(p), art::subRunFragment());

  return; 
}

DEFINE_ART_MODULE(sbn::SBNDBNBEXTRetriever)
