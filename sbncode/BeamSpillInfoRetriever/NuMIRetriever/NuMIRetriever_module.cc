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
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;
  
private:
  // input labels
  std::vector< sbn::NuMISpillInfo > fOutbeamInfos;
  double fTimePad;
  double fBFPEpsilion;
  std::string fURL;
  //MWRData mwrdata;
  std::string raw_data_label_;
  std::string fDeviceUsedForTiming;
  int TotalBeamSpills;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
};

sbn::NuMIRetriever::NuMIRetriever(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTimePad(p.get<double>("TimePadding", 0.5)), //epsilon in seconds, buffer time to look for spills
  fBFPEpsilion(p.get<double>("BFPEpsilon", 0.02)), // 20 ms, tuned for BNB, check for NuMI here might need to be larger
  raw_data_label_(p.get<std::string>("raw_data_label")),
  fDeviceUsedForTiming(p.get<std::string>("DeviceUsedForTiming")),
  bfp(ifbeam_handle->getBeamFolder(p.get<std::string>("Bundle"), p.get<std::string>("URL"), p.get<double>("TimeWindow")))
{

  bfp->set_epsilon(fBFPEpsilion);
  bfp->setValidWindow(500.);
  produces<std::vector<sbn::NuMISpillInfo>, art::InSubRun>();
  TotalBeamSpills = 0;
}

void sbn::NuMIRetriever::produce(art::Event &e)
{

  // If this is the first event in the run, then ignore it
  // We do not currently have the ability to figure out the first
  // spill that the DAQ was sensitive to, so don't try to save any
  // spill information
  //
  // TODO: long-term goal -- can we fix this?
  if (e.event() == 1) return;

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

    t_current_event = static_cast<double>(artdaq_ts)/(1000000000.); //check this offset... 
    if(gate_type == 2)
      t_previous_event = (static_cast<double>(frag.getLastTimestampNuMI()))/(1000000000.);
    else
      t_previous_event = (static_cast<double>(frag.getLastTimestampOther()))/(1000000000.);

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
    try{auto cur_vec_temp = bfp->GetNamedVector((t_previous_event)-fTimePad,"E:HP121[]");} catch (WebAPIException &we) {}
    
    try{auto cur_vec_temp_2 = bfp->GetNamedVector((t_current_event)+fTimePad,"E:VP121[]");} catch (WebAPIException &we) {}
    try{auto packed_MTGTDS_temp = bfp->GetNamedVector((t_current_event)+fTimePad, "E:MTGTDS[]");} catch(WebAPIException &we) {}
    std::vector<double> times_temps = bfp->GetTimeList(fDeviceUsedForTiming);

    int spill_count = 0;
    // Iterating through each of the beamline times
    for (size_t i = 0; i < times_temps.size(); i++) {

      // Only continue if these times are matched to our DAQ time
      // plus or minus some time padding, currently using 3.3 ms
      // which is half the Booster Rep Rate
      if(e.event() != 1){//We already addressed the "first event" above 
        if(times_temps[i] > t_current_event){continue;}
        if(times_temps[i] <= t_previous_event){continue;}
      }

      //count found spills
      spill_count++;

      //initialize all devices found in NuMISpillInfo.h in sbnobj
      double HRNDIR = -1.;
      double NSLINA = -1.;
      double NSLINB = -1.;
      double NSLINC = -1.;
      double NSLIND = -1.;
      double TOR101 = -1.;
      double TORTGT = -1.;
      double TR101D = -1.;
      double TRTGTD = -1.;
      std::vector< double > HP121;
      std::vector< double > VP121;
      std::vector< double > HPTGT;
      std::vector< double > VPTGT;
      std::vector< double > HITGT;
      std::vector< double > VITGT;
      std::vector< double > MTGTDS;
      double TRTGTD_time = -1.;
      std::cout << "Grabbing IFBeam info!" << std::endl;
      try{bfp->GetNamedData(times_temps[i], "E:TRTGTD@",&TRTGTD,&TRTGTD_time);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:TR101D",&TR101D);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:HRNDIR",&HRNDIR);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLINA",&NSLINA);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLINB",&NSLINB);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLINC",&NSLINC);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:NSLIND",&NSLIND);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:TOR101",&TOR101);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:TORTGT",&TORTGT);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      // BPM Positions and Intensities - they each have 7 elements
      // First is an average value of a few batches (often 2,3,4)
      // used for auto-tuning, so we should disregard it
      
      try{HP121 = bfp->GetNamedVector(times_temps[i], "E:HP121[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{VP121 = bfp->GetNamedVector(times_temps[i], "E:VP121[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      try{HPTGT = bfp->GetNamedVector(times_temps[i], "E:HPTGT[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      try{VPTGT = bfp->GetNamedVector(times_temps[i], "E:VPTGT[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      try{HITGT = bfp->GetNamedVector(times_temps[i], "E:HITGT[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      try{VITGT = bfp->GetNamedVector(times_temps[i], "E:VITGT[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      try{MTGTDS = bfp->GetNamedVector(times_temps[i], "E:MTGTDS[]");}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " <<"got exception: " << we.what() << "\n";}
      
      std::cout << "Finished getting IFBeam info" << std::endl;
      std::cout << "BFP Time: " << times_temps[i] << " TOROID Time: " << TRTGTD_time << " TOROID COUNT: " << TRTGTD << std::endl;
      unsigned long int time_closest_int = (int) TRTGTD_time;
      double time_closest_ns = (TRTGTD_time - time_closest_int)*1000000000;

      sbn::NuMISpillInfo NuMIbeamInfo;
      NuMIbeamInfo.TORTGT = TORTGT*1e12; //include factor of 1e12 protons in POT calculation
      NuMIbeamInfo.TOR101 = TOR101*1e12; //include factor of 1e12 protons in POT calculation
      NuMIbeamInfo.TRTGTD = TRTGTD;
      NuMIbeamInfo.TR101D = TR101D;
      NuMIbeamInfo.HRNDIR = HRNDIR;
      NuMIbeamInfo.NSLINA = NSLINA;
      NuMIbeamInfo.NSLINB = NSLINB;
      NuMIbeamInfo.NSLINC = NSLINC;
      NuMIbeamInfo.NSLIND = NSLIND;
      
      NuMIbeamInfo.HP121 = HP121;
      NuMIbeamInfo.VP121 = VP121;
      NuMIbeamInfo.HPTGT = HPTGT;
      NuMIbeamInfo.VPTGT = VPTGT;
      NuMIbeamInfo.HITGT = HITGT;
      NuMIbeamInfo.VITGT = VITGT;
      NuMIbeamInfo.MTGTDS = MTGTDS;

      NuMIbeamInfo.time = times_temps[i];
      NuMIbeamInfo.event = e.event();
      NuMIbeamInfo.spill_time_s = time_closest_int;
      NuMIbeamInfo.spill_time_ns = time_closest_ns;
      // Save the Number of DAQ Gates in the first saved spill
      if (spill_count == 1) {
        NuMIbeamInfo.daq_gates = number_of_gates_since_previous_event;
      }
      else {
        NuMIbeamInfo.daq_gates = 0;
      }
      
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

  // Clear the (now old) infos
  fOutbeamInfos.clear();

  return;
}

DEFINE_ART_MODULE(sbn::NuMIRetriever)


