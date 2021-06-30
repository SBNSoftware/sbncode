/*
 * NuMI Beam Spill Info Retriever module for ICARUS
 * Based heavily on code by Z. Pavlovic written for MicroBooNE
 * Based heavily on code by NOvA collaboration (Thanks NOvA!)
*/

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
#include "sbnobj/Common/POTAccounting/NuMISpillInfo.h"

#include "IFBeam_service.h"
#include "ifbeam_c.h"
//#include "MWRData.h"

#include <memory>
#include <optional>
#include <vector>

namespace sbn {
  class NuMIRetriever;
}

class sbn::NuMIRetriever : public art::EDProducer {
public:
  explicit NuMIRetriever(fhicl::ParameterSet const &p);
  NuMIRetriever(NuMIRetriever const&) = delete;
  NuMIRetriever(NuMIRetriever &&) = delete;
  NuMIRetriever& operator=(NuMIRetriever const&) = delete;
  NuMIRetriever& operator=(NuMIRetriever&&) = delete;
  
  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr);
  void endSubRun(art::SubRun& sr);
  
private:
  // input labels
  std::vector< sbn::NuMISpillInfo > fOutbeamInfos;
  int fTimePad;
  std::string fURL;
  //MWRData mwrdata;
  std::string raw_data_label_;
  std::string fDeviceUsedForTiming;
  int TotalBeamSpills;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;
  

};

sbn::NuMIRetriever::NuMIRetriever(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTimePad(p.get<double>("TimePadding", 0.5)), //epsilon in seconds, buffer time to look for spills
  raw_data_label_(p.get<std::string>("raw_data_label")),
  fDeviceUsedForTiming(p.get<std::string>("DeviceUsedForTiming")),
  bfp(ifbeam_handle->getBeamFolder(p.get<std::string>("Bundle"), p.get<std::string>("URL"), p.get<double>("TimeWindow"))),
  bfp_mwr( ifbeam_handle->getBeamFolder(p.get< std::string >("MultiWireBundle"), p.get< std::string >("URL"), p.get< double >("MWR_TimeWindow")))
{

  bfp->set_epsilon(0.02); //20 ms, tuned for BNB, check for NuMI here might need to be larger
  bfp_mwr->set_epsilon(0.5); //Not tuned, start with BNB value which is also not tuned
  bfp_mwr->setValidWindow(86400);
  produces<std::vector<sbn::NuMISpillInfo>, art::InSubRun>();
  TotalBeamSpills = 0;
}

void sbn::NuMIRetriever::produce(art::Event &e)
{
  int gate_type = 0;
  art::Handle< std::vector<artdaq::Fragment> > raw_data_ptr;
  e.getByLabel(raw_data_label_, "ICARUSTriggerUDP", raw_data_ptr);
  auto const & raw_data = (*raw_data_ptr);

  double t_current_event  = 0;
  double t_previous_event = 0;
  double number_of_gates_since_previous_event = 0;

  for(auto raw_datum : raw_data){

    uint64_t artdaq_ts = raw_datum.timestamp();
    icarus::ICARUSTriggerUDPFragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);
    gate_type = datastream_info.gate_type;
    number_of_gates_since_previous_event = frag.getDeltaGatesNuMI();

    t_current_event = static_cast<double>(artdaq_ts)/(1000000000); //check this offset... 
    if(gate_type == 2)
      t_previous_event = (static_cast<double>(frag.getLastTimestampNuMI()))/(1000000000);
    else
      t_previous_event = (static_cast<double>(frag.getLastTimestampOther()))/(1000000000);

  }

  std::cout << std::setprecision(19) << "Previous : " << t_previous_event << ", Current : " << t_current_event << std::endl;
  //We only want to process NuMI gates, i.e. type 2
  if(gate_type == 2)
  {
    // Keep track of the number of beam gates the DAQ thinks
    //   are in this job   
    TotalBeamSpills += number_of_gates_since_previous_event;
    
    // These lines get everything primed within the IFBeamDB
    //   They seem redundant but they are needed
    std::cout << "Trying beam info" << std::endl;
    try{auto cur_vec_temp = bfp->GetNamedVector((t_previous_event)-fTimePad,"E:HP121[]");} catch (WebAPIException &we) {}
    
    //try{auto cur_vec_temp_2 = bfp->GetNamedVector((t_current_event)+fTimePad,"E:VP121[]");} catch (WebAPIException &we) {}
    //try{auto packed_MTGTDS_temp = bfp->GetNamedVector((t_current_event)+fTimePad, "E:MTGTDS[]");} catch(WebAPIException &we) {}
    std::cout << "IFBeam test succeeded" << std::endl;
    std::vector<double> times_temps = bfp->GetTimeList(fDeviceUsedForTiming);

    int spill_count = 0;
    // Iterating through each of the beamline times
    for (size_t i = 0; i < times_temps.size(); i++) {

      // Only continue if these times are matched to our DAQ time
      // plus or minus some time padding, currently using 3.3 ms
      // which is half the Booster Rep Rate
      if(e.event() != 1){//We already addressed the "first event" above 
        if(times_temps[i] > (t_current_event+fTimePad)){continue;}
        if(times_temps[i] <= (t_previous_event-fTimePad)){continue;}
      }

      //count found spills
      spill_count++;

      //initialize all devices found in NuMISpillInfo.h in sbnobj
      float HRNDIR;
      float NSLINA;
      float NSLINB;
      float NSLINC;
      float NSLIND;
      float TR101D;
      float TRTGTD;
      std::vector< double > HP121;
      std::vector< double > VP121;
      std::vector< double > HPTGT;
      std::vector< double > VPTGT;
      std::vector< double > HITGT;
      std::vector< double > VITGT;
      std::vector< double > MTGTDS;
      float TRTGTD_time;
      std::cout << "Grabbing IFBeam info!" << std::endl;
      try{bfp->GetNamedData(times_temps[i], "E:TRTGTD@",&TRTGTD,&TRTGTD_time);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:TR101D",&TR101D);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:HRNDIR",&HRNDIR);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLINA",&NSLINA);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLINB",&NSLINB);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLINC",&NSLINC);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLIND",&NSLIND);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      // BPM Positions and Intensities - they each have 7 elements
      // First is an average value of a few batches (often 2,3,4)
      // used for auto-tuning, so we should disregard it
      
      //try{HP121 = bfp->GetNamedVector(times_temps[i], "E:HP121[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      //try{VP121 = bfp->GetNamedVector(times_temps[i], "E:VP121[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      //try{HPTGT = bfp->GetNamedVector(times_temps[i], "E:HPTGT");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      //try{VPTGT = bfp->GetNamedVector(times_temps[i], "E:VPTGT");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      //try{HITGT = bfp->GetNamedVector(times_temps[i], "E:HITGT");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      //try{VITGT = bfp->GetNamedVector(times_temps[i], "E:VITGT");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      
      
      HP121 = bfp->GetNamedVector(times_temps[i], "E:HP121[]");
      std::cout << "Got first device!" << std::endl;
      //VP121 = bfp->GetNamedVector(times_temps[i], "E:VP121[]");
      //std::cout << "Got second device!" << std::endl;
      //HPTGT = bfp->GetNamedVector(times_temps[i], "E:HPTGT[]");
      //std::cout << "Got third device!" << std::endl;
      //VPTGT = bfp->GetNamedVector(times_temps[i], "E:VPTGT[]");
      //HITGT = bfp->GetNamedVector(times_temps[i], "E:HITGT[]");
      //VITGT = bfp->GetNamedVector(times_temps[i], "E:VITGT[]");
      
      std::cout << "Finished getting IFBeam info" << std::endl;
      //unsigned long int time_closest_int = (int) TRTGTD_time;
      //double time_closest_ns = (TRTGTD_time - time_closest_int)*1000000000;
      sbn::NuMISpillInfo NuMIbeamInfo;
      NuMIbeamInfo.TRTGTD = TRTGTD;
      NuMIbeamInfo.TR101D = TR101D;
      NuMIbeamInfo.HRNDIR = HRNDIR;
      NuMIbeamInfo.NSLINA = NSLINA;
      NuMIbeamInfo.NSLINB = NSLINB;
      NuMIbeamInfo.NSLINC = NSLINC;
      NuMIbeamInfo.NSLIND = NSLIND;
      
      NuMIbeamInfo.HP121 = HP121;
      //NuMIbeamInfo.VP121 = VP121;
      //NuMIbeamInfo.HPTGT = HPTGT;
      //NuMIbeamInfo.VPTGT = VPTGT;
      //NuMIbeamInfo.HITGT = HITGT;
      //NuMIbeamInfo.VITGT = VITGT;
      
      fOutbeamInfos.push_back(NuMIbeamInfo);
    }
    if(spill_count > number_of_gates_since_previous_event)
      std::cout << "Event Spills : " << spill_count << ", DAQ Spills : " << number_of_gates_since_previous_event << " \t \t ::: WRONG!"<< std::endl;
    else
      std::cout << "Event Spills : " << spill_count << ", DAQ Spills : " << number_of_gates_since_previous_event << std::endl;
  }
}

void sbn::NuMIRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

void sbn::NuMIRetriever::endSubRun(art::SubRun& sr)
{
  // We will add all of the BNBSpillInfo data-products to the
  // art::SubRun so it persists
  // currently this is ~2.7 kB/event or ~0.07 kB/spill

  std::cout << "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
  std::cout << "Total number of Selected Spills : " << fOutbeamInfos.size() << std::endl;

  auto p =  std::make_unique< std::vector< sbn::NuMISpillInfo > >(fOutbeamInfos);

  sr.put(std::move(p), art::subRunFragment());

  return;
}

DEFINE_ART_MODULE(sbn::NuMIRetriever)


