////////////////////////////////////////////////////////////////////////
// Class:       EXTRetriever
// Plugin Type: producer 
// File:        EXTRetriever_module.cc
//
// Created by hand Thurs June 24th 2021 by J. Zennamo (FNAL)
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
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcorealg/Geometry/Exceptions.h"

#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerUDPFragment.hh"

#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"

#include <memory>
#include <optional>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <time.h>

namespace sbn {
  class EXTRetriever;
}

class sbn::EXTRetriever : public art::EDProducer {
public:
  explicit EXTRetriever(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EXTRetriever(EXTRetriever const&) = delete;
  EXTRetriever(EXTRetriever&&) = delete;
  EXTRetriever& operator=(EXTRetriever const&) = delete;
  EXTRetriever& operator=(EXTRetriever&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

private:
  std::vector< sbn::EXTCountInfo > fOutExtInfos;
 
  // input labels
  std::string raw_data_label_;
  int TotalEXTCounts;  

};


sbn::EXTRetriever::EXTRetriever(fhicl::ParameterSet const& p)
  : EDProducer{p},
  raw_data_label_(p.get<std::string>("raw_data_label"))
{
 
  produces< std::vector< sbn::EXTCountInfo >, art::InSubRun >();
  TotalEXTCounts = 0;
}

void sbn::EXTRetriever::produce(art::Event& e)
{
  
  //Here we read in the artdaq Fragments and extract three pieces of information:
  // 1. The time of the current event, t_current_event
  // 2. the time of the previously triggered event, t_previous_event (NOTE: Events are non-sequential!)
  // 3. the number of beam spills since the previously triggered event, number_of_gates_since_previous_event
  
  int gate_type = 0;
  art::Handle< std::vector<artdaq::Fragment> > raw_data_ptr;
  e.getByLabel(raw_data_label_, "ICARUSTriggerUDP", raw_data_ptr);
  auto const & raw_data = (*raw_data_ptr);

  double number_of_gates_since_previous_event = 0;
  
  for(auto raw_datum : raw_data){
   
    icarus::ICARUSTriggerUDPFragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);
    gate_type = datastream_info.gate_type;
    number_of_gates_since_previous_event = frag.getDeltaGatesBNB();
  
  }
  
  //We only want to process EXT gates, i.e. type 3
  if(gate_type == 3)
  {
    // Keep track of the number of beam gates the DAQ thinks 
    //   are in this file
    TotalEXTCounts += number_of_gates_since_previous_event;
   
      //Store everything in our data-product
      sbn::EXTCountInfo extInfo;
      extInfo.gates_since_last_trigger = number_of_gates_since_previous_event;

      fOutExtInfos.push_back(extInfo);
      // We do not write these to the art::Events because 
      // we can filter events but want to keep all the POT 
      // information, so we'll write it to the SubRun
    
  }//end check on gate type
} //end loop over events


void sbn::EXTRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

//____________________________________________________________________________                                                                                                                                                                                      
void sbn::EXTRetriever::endSubRun(art::SubRun& sr)
{
  // We will add all of the EXTCountInfo data-products to the 
  // art::SubRun so it persists 
  auto p =  std::make_unique< std::vector< sbn::EXTCountInfo > >(fOutExtInfos);

  sr.put(std::move(p));

  return;
}

DEFINE_ART_MODULE(sbn::EXTRetriever)    
