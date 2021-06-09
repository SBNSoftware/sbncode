////////////////////////////////////////////////////////////////////////
// Class:       Retriever
// Plugin Type: producer 
// File:        Retriever_module.cc
//
// Created by hand Wed April 9 2021 by J. Zennamo (FNAL)
// Based heavily on code by Z. Pavlovic written for MicroBooNE 
// Based heavily on code by NOvA collaboration (Thanks NOvA!):
//       https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/BNBInfo_module.cc
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
#include "sbnobj/Common/POTAccounting/BeamSpillInfo.h"

#include "IFBeam_service.h"
#include "ifbeam_c.h"
#include "MWRData.h"

#include <memory>
#include <optional>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <time.h>

namespace sbn {
  class Retriever;
}

class sbn::Retriever : public art::EDProducer {
public:
  explicit Retriever(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Retriever(Retriever const&) = delete;
  Retriever(Retriever&&) = delete;
  Retriever& operator=(Retriever const&) = delete;
  Retriever& operator=(Retriever&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr);
  void endSubRun(art::SubRun& sr);

private:
  // input labels
  std::vector< sbn::BeamSpillInfo > fOutbeamInfos;
  int fTimePad;
  std::string fURL;
  MWRData mwrdata;
  std::string raw_data_label_;
 
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;

};


sbn::Retriever::Retriever(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTimePad(p.get<double>("TimePadding",0.0333)), //seconds 
  raw_data_label_(p.get<std::string>("raw_data_label")),
  bfp(     ifbeam_handle->getBeamFolder(p.get< std::string >("Bundle"), p.get< std::string >("URL"), p.get< double >("TimeWindow"))),
  bfp_mwr( ifbeam_handle->getBeamFolder(p.get< std::string >("MultiWireBundle"), p.get< std::string >("URL"), p.get< double >("TimeWindow")))
{
  bfp->set_epsilon(0.0333); // how close in time does the spill time have to be from the DAQ time (in seconds).
  bfp_mwr->set_epsilon(100); // how close in time does the spill time have to be from the DAQ time (in seconds).
  produces< std::vector< sbn::BeamSpillInfo >, art::InSubRun >();
  
}

void sbn::Retriever::produce(art::Event& e)
{

  //Get the time of the event and the previous event
  // time since epoch sec . nanoseconds
  //1616202843.932767280
  //.932767280
  //.030000000
  
  int gate_type = 0;
  art::Handle< std::vector<artdaq::Fragment> > raw_data_ptr;
  e.getByLabel(raw_data_label_, "ICARUSTriggerUDP", raw_data_ptr);
  auto const & raw_data = (*raw_data_ptr);
  //  double t_previous_event = 1616202844.932767280;
  //  double t_current_event  = 1616202843.932767280;

  double t_previous_event = 0;
  double t_current_event  = 0;
  std::cout << "This many raw_datums : " << raw_data.size() << std::endl;
  for(auto raw_datum : raw_data){
   
    uint64_t artdaq_ts = raw_datum.timestamp();
    icarus::ICARUSTriggerUDPFragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);  
    gate_type = datastream_info.gate_type;
    double delta_gates = frag.getDeltaGatesBNB();
    std::cout << "Gate_type: " << gate_type << std::endl;
    if(gate_type == 1)
      std::cout << "Number of beam gates since the last spill of this type: " << delta_gates << std::endl;
    std::cout << "Fragment timestamp: " << artdaq_ts << std::endl; 
    t_current_event = static_cast<double>(artdaq_ts)/(1000000000); //check this offset...
    if(gate_type == 1)
      t_previous_event = (static_cast<double>(frag.getLastTimestampBNB()))/(1e9);
    else
      t_previous_event = (static_cast<double>(frag.getLastTimestampOther()))/(1000000000);
    int gate_type = datastream_info.gate_type;
    if(gate_type != 1) continue;

    t_current_event = datastream_info.getNanoseconds_since_UTC_epoch() + 2; //check this offset...
    t_previous_event = frag.getLastTimestampBNB();
     
  }
  
  std::cout << "Previous : " << t_previous_event << ", Current : " << t_current_event << std::endl;
  
  //if(t_previous_event == t_current_event) return -1;
  
  try{
    double test_curr; double test_t1;
    bfp->GetNamedData((t_previous_event)-fTimePad,"E:THCURR",&test_curr,&test_t1);
    std::cout << "got values " << test_curr <<  "for E:THCURR at time " << test_t1 << "\n";
  }
  catch (WebAPIException &we) {
    std::cout << "got exception: " << we.what() << "\n";
  }
  std::cout.precision(17);


  try{
    auto packed_M876BB_temp = bfp_mwr->GetNamedVector((t_previous_event)-35,"E:M875BB{4440:888}.RAW");
  }
  catch (WebAPIException &we) {
    //we have to have this or it doesn't work...
>>>>>>> b831f5685458e11a49d3ee5445bd0fc65e47295c
  }
  
  std::cout << std::setprecision(19) << "Previous : " << t_previous_event << ", Current : " << t_current_event << std::endl;
  
  //if(t_previous_event == t_current_event) continue;

  if(gate_type == 1)
  {
    std::cout << "Is BNB event! Extracting POT!" << std::endl;
    try{
      double test_curr; double test_t1;
      bfp->GetNamedData((t_previous_event)-fTimePad,"E:THCURR",&test_curr,&test_t1);
      std::cout << "got values " << test_curr <<  "for E:THCURR at time " << test_t1 << "\n";
    }
    catch (WebAPIException &we) {
      std::cout << "got exception: " << we.what() << "\n";
    }
    std::cout.precision(17);
    
    
    try{
      auto packed_M876BB_temp = bfp_mwr->GetNamedVector((t_previous_event)-35,"E:M875BB{4440:888}.RAW");
    }
    catch (WebAPIException &we) {
      //we have to have this or it doesn't work...
    }
    
    //  bfp_mwr->setValidWindow(86400);  
    //std::cout << "Is this a valid Window " << bfp_mwr->getValidWindow()<<std::endl;
    
    
  //  std::vector< std::vector<std::string> > unpacked_M876BB_str;
  // unpacked_M876BB_str.resize(3);
    std::vector<std::string> unpacked_M876BB_str;
    std::string packed_data_str; 
    //  int dev = 0;
    std::vector<std::string> vars = bfp_mwr->GetDeviceList();
    double t_mwr;
    std::cout << int(vars.size()) << " MWR Devices " << std::endl;
    for (int i = 0; i < int(vars.size()); i++) {
      if(vars[i].empty()) continue;
      
      /// Check the device name and interate the double-vector index
      // if(vars[i].M875BB) dev = 0;
      //M876BB;
      //MMBTBB;
      //std::cout << vars[i] << std::endl;
      
      
    //var - > E:M875BB{888:888}.RAW
      t_mwr = 0;
      std::vector<double> packed_M876BB = bfp_mwr->GetNamedVector((t_previous_event)-35,vars[i],&t_mwr);
      
      packed_data_str.clear();
      packed_data_str += std::to_string(int(t_mwr));
      packed_data_str.append(",");
      packed_data_str.append(vars[i]);
      packed_data_str.append(",,");
      
      for(int j = 0; j < int(packed_M876BB.size()); j++){
	packed_data_str += std::to_string(int(packed_M876BB[j]));
	if(j < int(packed_M876BB.size())-1)
	  packed_data_str.append(",");
      }
      
      //    std::cout << packed_data_str << std::endl;
      
      auto unpacked_M876BB_str_temp = mwrdata.unpackMWR(packed_data_str,long(-35));
      unpacked_M876BB_str.insert(unpacked_M876BB_str.end(),
				 unpacked_M876BB_str_temp.begin(),
				 unpacked_M876BB_str_temp.end());
      
      
      //std::cout <<  int(unpacked_M876BB_str.size()) << " unpacked datas ... " << std::endl;
    }

    //for(int j = 0; j < int(unpacked_M876BB_str.size()); j++){
    //  std::cout << unpacked_M876BB_str[j] << std::endl;
  //}
    
    
    std::vector<double> times_temps = bfp->GetTimeList();
    for (size_t i = 0; i < times_temps.size(); i++) {
      if(times_temps[i] > (t_current_event+fTimePad)){break;}
      
      double TOR860 = 0;
      double TOR875 = 0;
      double LM875A = 0;
      double LM875B = 0;
      double LM875C = 0;
      double HP875 = 0;
      double VP875 = 0;
      double HPTG1 = 0;
      double VPTG1 = 0;
      double HPTG2 = 0;
      double VPTG2 = 0;
      double BTJT2 = 0;
      double THCURR = 0;
      
      double TOR860_time = 0;
      
      bfp->GetNamedData(times_temps[i], "E:TOR860@,E:TOR875,E:LM875A,E:LM875B,E:LM875C,E:HP875,E:VP875,E:HPTG1,E:VPTG1,E:HPTG2,E:VPTG2,E:BTJT2,E:THCURR",&TOR860,&TOR860_time,&TOR875,&LM875A,&LM875B,&LM875C,&HP875,&VP875,&HPTG1,&VPTG1,&HPTG2,&VPTG2,&BTJT2,&THCURR);
      
      
      unsigned long int time_closest_int = (int) TOR860_time;
      double time_closest_ns = (TOR860_time - time_closest_int)*1e9;
      
      sbn::BeamSpillInfo beamInfo;
      beamInfo.TOR860 = TOR860;
      beamInfo.TOR875 = TOR875;
      beamInfo.LM875A = LM875A;
      beamInfo.LM875B = LM875B;
      beamInfo.LM875C = LM875C;
      beamInfo.HP875 = HP875;
      beamInfo.VP875 = VP875;
      beamInfo.HPTG1 = HPTG1;
      beamInfo.VPTG1 = VPTG1;
      beamInfo.HPTG2 = HPTG2;
      beamInfo.VPTG2 = VPTG2;
      beamInfo.BTJT2 = BTJT2;
      beamInfo.THCURR = THCURR;
      beamInfo.spill_time_s = time_closest_int;
      beamInfo.spill_time_ns = time_closest_ns;    
      
      fOutbeamInfos.push_back(beamInfo);
    
    
    }
  }
}

void sbn::Retriever::beginSubRun(art::SubRun& sr)
{
  return;
}

//____________________________________________________________________________                                                                                                                                                                                      
void sbn::Retriever::endSubRun(art::SubRun& sr)
{

  auto p =  std::make_unique< std::vector< sbn::BeamSpillInfo > >(fOutbeamInfos);

  sr.put(std::move(p));

  return;
}

DEFINE_ART_MODULE(sbn::Retriever)    
