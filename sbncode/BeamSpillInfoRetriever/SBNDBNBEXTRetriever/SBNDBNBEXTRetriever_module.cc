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
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"
#include "sbncode/BeamSpillInfoRetriever/SBNDPOTTools.h"

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
  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  // input labels
  float TotalEXTCounts;  
};

sbn::SBNDBNBEXTRetriever::SBNDBNBEXTRetriever(fhicl::ParameterSet const & params)
  : EDProducer{params} {
  produces< std::vector< sbn::EXTCountInfo >, art::InSubRun >();
  TotalEXTCounts = 0;
}

void sbn::SBNDBNBEXTRetriever::produce(art::Event & e)
{

  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);
  TotalEXTCounts += triggerInfo.number_of_gates_since_previous_event;
  //Store everything in our data-product
  sbn::EXTCountInfo extInfo;
  extInfo.gates_since_last_trigger = triggerInfo.number_of_gates_since_previous_event;
  fOutExtInfos.push_back(extInfo);
}

sbn::TriggerInfo_t sbn::SBNDBNBEXTRetriever::extractTriggerInfo(art::Event const& e) const {
  // Using TDC for current event, but PTB for previous event.
  // Don't worry about BESOffset since we only record number of gates
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);

  art::InputTag TDC_itag("daq", "ContainerTDCTIMESTAMP");
  auto TDC_cont_frags = e.getHandle<artdaq::Fragments>(TDC_itag);

  PTBInfo_t PTBInfo;
  TriggerInfo_t triggerInfo;
  PTBInfo = extractPTBInfo(PTB_cont_frags, 4);

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

  double PTBandCurrOffset = PTBInfo.currPTBTimeStamp - triggerInfo.t_current_event;

  // Catch for an issue seen a few times where PTB off by a second.
  // Only need to correct prevTS because either currTS is from TDC
  // or there is no offset between currTS and PTB.
  if(abs(PTBandCurrOffset) >= 0.5){
    triggerInfo.t_previous_event-=PTBandCurrOffset;
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << "Offset between PTB and TDC, " << PTBandCurrOffset << std::endl;
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << "Corrected previous event TS is " << std::setprecision(19) << triggerInfo.t_previous_event << std::endl;
  }

  return triggerInfo;
}

void sbn::SBNDBNBEXTRetriever::beginSubRun(art::SubRun& sr)
{
  TotalEXTCounts = 0;
  fOutExtInfos = {};
  return;
}

void sbn::SBNDBNBEXTRetriever::endSubRun(art::SubRun& sr)
{
  // We will add all of the EXTCountInfo data-products to the 
  // art::SubRun so it persists 

  mf::LogDebug("SBNDBNBEXTRetriever")<< "Total number of DAQ Spills : " << TotalEXTCounts << std::endl;
  auto p =  std::make_unique< std::vector< sbn::EXTCountInfo > >(fOutExtInfos);

  sr.put(std::move(p), art::subRunFragment());

  return; 
}

DEFINE_ART_MODULE(sbn::SBNDBNBEXTRetriever)
