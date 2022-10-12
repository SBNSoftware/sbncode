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
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcorealg/Geometry/Exceptions.h"

#include "artdaq-core/Data/Fragment.hh"
#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerV2Fragment.hh"

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
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<std::string> RawDataLabel {
      Name{ "raw_data_label" },
      Comment{ "art data product instance name for trigger information (product label is 'daq')" }
      };
    
  }; // Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  explicit EXTRetriever(Parameters const& params);
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


sbn::EXTRetriever::EXTRetriever(Parameters const& params)
  : EDProducer{params},
  raw_data_label_(params().RawDataLabel())
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
  int trigger_type = 0;
  auto const & raw_data = e.getProduct< std::vector<artdaq::Fragment> >({ raw_data_label_, "ICARUSTriggerV2" });

  unsigned int number_of_gates_since_previous_event = 0;
  bool isBNBOffBeam = false;
  bool isNuMIOffBeam = false;
  bool isMajority = false;
  bool isMinBias = false;

  for(auto raw_datum : raw_data){
   
    icarus::ICARUSTriggerV2Fragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerV2String(buffer);
    gate_type = datastream_info.gate_type;
    trigger_type = datastream_info.trigger_type;
    
    if(gate_type == 3)
    {
      isBNBOffBeam = true;
      number_of_gates_since_previous_event = frag.getDeltaGatesBNBOff();
    }
    else if(gate_type == 4)
    {
      isNuMIOffBeam = true;
      number_of_gates_since_previous_event = frag.getDeltaGatesNuMIOff();
    }
    else 
    {
      throw art::Exception(art::errors::StdException) << "Unsupported gate type for EXTRetriever Module: " << gate_type << "! Aborting";
    }
    if(trigger_type == 0)
      isMajority = true;
    if(trigger_type == 1)
      isMinBias = true;
  }
  
  //We only want to process EXT gates, i.e. type 3 (BNBOffBeam) or 4 (NuMIOffBeam) 
  if(gate_type == 3 || gate_type == 4)
  {
    // Keep track of the number of beam gates the DAQ thinks 
    //   are in this file
    TotalEXTCounts += number_of_gates_since_previous_event;
   
      //Store everything in our data-product
      sbn::EXTCountInfo extInfo;
      extInfo.gates_since_last_trigger = number_of_gates_since_previous_event;
      extInfo.isBNBOffBeam = isBNBOffBeam;
      extInfo.isNuMIOffBeam = isNuMIOffBeam;
      extInfo.isMajority = isMajority;
      extInfo.isMinBias = isMinBias;
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
