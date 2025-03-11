////////////////////////////////////////////////////////////////////////
// Class:       SBNDBNBEXTRetriever
// Plugin Type: producer 
// File:        SBNDBNBEXTRetriever_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"
#include "sbncode/BeamSpillInfoRetriever/SBNDPOTTools.h"

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
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);
  PTBInfo_t PTBInfo = extractPTBInfo(PTB_cont_frags, 4);
  int SingleEventGateCounter = PTBInfo.GateCounter;
 
  TotalEXTCounts += SingleEventGateCounter;
  sbn::EXTCountInfo extInfo;
  extInfo.gates_since_last_trigger = SingleEventGateCounter;
  fOutExtInfos.push_back(extInfo);
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
