////////////////////////////////////////////////////////////////////////
// Class:       BNBEXTRetriever
// Plugin Type: producer 
// File:        BNBEXTRetriever_module.cc
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
#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerV3Fragment.hh"

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
  class BNBEXTRetriever;
}

class sbn::BNBEXTRetriever : public art::EDProducer {
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
  
  
  explicit BNBEXTRetriever(Parameters const& params);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BNBEXTRetriever(BNBEXTRetriever const&) = delete;
  BNBEXTRetriever(BNBEXTRetriever&&) = delete;
  BNBEXTRetriever& operator=(BNBEXTRetriever const&) = delete;
  BNBEXTRetriever& operator=(BNBEXTRetriever&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

private:
  std::vector< sbn::EXTCountInfo > fOutExtInfos;
 
  // input labels
  std::string raw_data_label_;
  float TotalEXTCounts;  
  float scale_factor;

};


sbn::BNBEXTRetriever::BNBEXTRetriever(Parameters const& params)
  : EDProducer{params},
  raw_data_label_(params().RawDataLabel())
{
 
  produces< std::vector< sbn::EXTCountInfo >, art::InSubRun >();
  TotalEXTCounts = 0;
  scale_factor = 0;
}

void sbn::BNBEXTRetriever::produce(art::Event& e)
{
  
  //Here we read in the artdaq Fragments and extract three pieces of information:
  // 1. The time of the current event, t_current_event
  // 2. the time of the previously triggered event, t_previous_event (NOTE: Events are non-sequential!)
  // 3. the number of beam spills since the previously triggered event, number_of_gates_since_previous_event
  
  int gate_type = 0;
  auto const & raw_data = e.getProduct< std::vector<artdaq::Fragment> >({ raw_data_label_, "ICARUSTriggerV3" });

  unsigned int number_of_gates_since_previous_event = 0;
  
  for(auto raw_datum : raw_data){
   
    icarus::ICARUSTriggerV3Fragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerV3String(buffer);
    gate_type = datastream_info.gate_type;
    number_of_gates_since_previous_event = frag.getDeltaGatesBNBOffMaj();
    scale_factor = 1. - (1./frag.getDeltaGatesBNBOffMinbias());
    std::cout << "BNB OFF MAJ : " << frag.getDeltaGatesBNBOffMaj() << std::endl; 
    std::cout << "Scale Factor : " << scale_factor << std::endl; 

  }
  
  //We only want to process EXT gates, i.e. type 3
  if(gate_type == 3)
  {
    // Keep track of the number of beam gates the DAQ thinks 
    //   are in this file
    TotalEXTCounts += number_of_gates_since_previous_event*scale_factor;
   
      //Store everything in our data-product
      sbn::EXTCountInfo extInfo;
      extInfo.gates_since_last_trigger = number_of_gates_since_previous_event;

      fOutExtInfos.push_back(extInfo);
      // We do not write these to the art::Events because 
      // we can filter events but want to keep all the POT 
      // information, so we'll write it to the SubRun
    
  }//end check on gate type
} //end loop over events


void sbn::BNBEXTRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

//____________________________________________________________________________                                                                                                                                                                                      
void sbn::BNBEXTRetriever::endSubRun(art::SubRun& sr)
{
  // We will add all of the EXTCountInfo data-products to the 
  // art::SubRun so it persists 
  auto p =  std::make_unique< std::vector< sbn::EXTCountInfo > >(fOutExtInfos);

  sr.put(std::move(p), art::subRunFragment());

  return;
}

DEFINE_ART_MODULE(sbn::BNBEXTRetriever)    
