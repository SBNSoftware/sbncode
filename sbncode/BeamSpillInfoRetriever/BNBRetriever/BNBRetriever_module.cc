/** ********************************************************************
 * @file BNBRetriever_module.cc
 * @date Wed April 9 2021
 * @author J. Zennamo (FNAL)
 * 
 * Based heavily on code by Z. Pavlovic written for MicroBooNE 
 * Based heavily on code by NOvA collaboration (Thanks NOvA!):
 *        https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/BNBInfo_module.cc
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
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"

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
  class BNBRetriever;
}

class sbn::BNBRetriever : public art::EDProducer {
public:
  explicit BNBRetriever(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BNBRetriever(BNBRetriever const&) = delete;
  BNBRetriever(BNBRetriever&&) = delete;
  BNBRetriever& operator=(BNBRetriever const&) = delete;
  BNBRetriever& operator=(BNBRetriever&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr);
  void endSubRun(art::SubRun& sr);

private:
  // input labels
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  int fTimePad;
  std::string fURL;
  MWRData mwrdata;
  std::string raw_data_label_;
  std::string fDeviceUsedForTiming;
  int TotalBeamSpills;  
  //
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;



};


sbn::BNBRetriever::BNBRetriever(fhicl::ParameterSet const& p)
  : EDProducer{p},
  fTimePad(p.get<double>("TimePadding",0.0333)), //seconds 
  raw_data_label_(p.get<std::string>("raw_data_label")),
  fDeviceUsedForTiming(p.get<std::string>("DeviceUsedForTiming")),
  bfp(     ifbeam_handle->getBeamFolder(p.get< std::string >("Bundle"), p.get< std::string >("URL"), p.get< double >("TimeWindow"))),
  bfp_mwr( ifbeam_handle->getBeamFolder(p.get< std::string >("MultiWireBundle"), p.get< std::string >("URL"), p.get< double >("MWR_TimeWindow")))
{
 
  // how close in time does the spill time have to be from the DAQ time (in seconds).
  // If these are too large then it fails to capture the device 
  // If these are too small then the time jitter in devices means we miss good data 
  bfp->set_epsilon(0.02); //20 ms, this was tuned by hand and compared to IFBeamDB times  
  
  bfp_mwr->set_epsilon(0.5); // TO BE TUNED!
  bfp_mwr->setValidWindow(86400);  
  produces< std::vector< sbn::BNBSpillInfo >, art::InSubRun >();
  TotalBeamSpills = 0;
}

void sbn::BNBRetriever::produce(art::Event& e)
{
  
  //Here we read in the artdaq Fragments and extract three pieces of information:
  // 1. The time of the current event, t_current_event
  // 2. the time of the previously triggered event, t_previous_event (NOTE: Events are non-sequential!)
  // 3. the number of beam spills since the previously triggered event, number_of_gates_since_previous_event
  
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
    number_of_gates_since_previous_event = frag.getDeltaGatesBNB();
  
    t_current_event = static_cast<double>(artdaq_ts)/(1000000000); //check this offset...
    if(gate_type == 1)
      t_previous_event = (static_cast<double>(frag.getLastTimestampBNB()))/(1e9);
    else
      t_previous_event = (static_cast<double>(frag.getLastTimestampOther()))/(1000000000);
    
  }
  
  std::cout << std::setprecision(19) << "Previous : " << t_previous_event << ", Current : " << t_current_event << std::endl;

  //We only want to process BNB gates, i.e. type 1 
  if(gate_type == 1)
  {
    // Keep track of the number of beam gates the DAQ thinks 
    //   are in this file
    TotalBeamSpills += number_of_gates_since_previous_event;
   
    // These lines get everything primed within the IFBeamDB
    //   They seem redundant but they are needed
    try{auto cur_vec_temp = bfp->GetNamedVector((t_previous_event)-fTimePad,"E:THCURR");} catch (WebAPIException &we) {}      
    try{auto packed_M876BB_temp = bfp_mwr->GetNamedVector((t_current_event)+fTimePad,"E:M875BB{4440:888}.RAW");} catch (WebAPIException &we) {}


    //The multiwire chambers provide their
    // data in a vector format but we'll have 
    // to sort through it in std::string format
    // to correctly unpack it
    std::vector< std::vector< std::vector< int > > >  unpacked_MWR;
    std::vector< std::vector< double> > MWR_times;
    unpacked_MWR.resize(3);
    MWR_times.resize(3);
    std::string packed_data_str; 
    
    //Create a list of all the MWR devices with their different
    // memory buffer increments 
    // generally in the format: "E:<Device>.{Memory Block}"
    std::vector<std::string> vars = bfp_mwr->GetDeviceList();
    
    // Tracking the time from the IFBeamDB
    double time_for_mwr;    
    
    // this is an iterator to track which of the 
    // three devices we will be working with
    int dev = 0;
    
    // The MWR devices are annoying and have confusing buffer
    // what we'll do is sort through all of them first and then 
    // match them to the closest spills in time
    // 
    int t_steps = int(fabs((t_previous_event - fTimePad) - (t_current_event + fTimePad))/0.5)+25;
    
    for(int t = 0; t < t_steps; t++){//Iterate through time increments
      for (int i = 0; i < int(vars.size()); i++) {// Iterate through the devices

	//Make sure we have a device
	if(vars[i].empty()) continue;
	
	/// Check the device name and interate the double-vector index
	if(vars[i].find("M875BB") != std::string::npos ) dev = 0;
	else if(vars[i].find("M876BB") != std::string::npos ) dev = 1;
	else if(vars[i].find("MMBTBB") != std::string::npos ) dev = 2;
	else{continue;}

	time_for_mwr = 0;
	
	try{
	  //Pull the MWR data for the device
	  // these data are "packed"
	  std::vector<double> packed_MWR = bfp_mwr->GetNamedVector((t_previous_event)-fTimePad+double(0.5*t),vars[i],&time_for_mwr);

	  //We'll convert this into a format
	  // that we can unpack doubles >> strings
	  //
	  packed_data_str.clear();
	  packed_data_str += std::to_string(int(time_for_mwr));
	  packed_data_str.append(",");
	  packed_data_str.append(vars[i]);
	  packed_data_str.append(",,");
	  
	  for(int j = 0; j < int(packed_MWR.size()); j++){
	    packed_data_str += std::to_string(int(packed_MWR[j]));
	    if(j < int(packed_MWR.size())-1)
	      packed_data_str.append(",");
	  }
	  
	  // Use Zarko's unpacking function to turn this into consumeable data
	  std::vector<double> MWR_times_temp;

	  // There is a 35 ms offset between the toriod and the MWR times
	  //   we'll just remove that here to match to the spill times
	  std::vector< std::vector< int > > unpacked_MWR_temp = mwrdata.unpackMWR(packed_data_str,MWR_times_temp,-0.035);
	  
	  //There are four events that are packed into one MWR IFBeam entry
	  for(int s = 0; s < int(unpacked_MWR_temp.size()); s++){
	  
	    // If this entry has a unique time them store it for later	  
	    if(std::find(MWR_times[dev].begin(), MWR_times[dev].end(), MWR_times_temp[s]) == MWR_times[dev].end()){
	      unpacked_MWR[dev].push_back(unpacked_MWR_temp[s]);
	      MWR_times[dev].push_back(MWR_times_temp[s]);
	    }//check for unique time 

	  }//Iterate through the unpacked events
	}//try
	catch (WebAPIException &we) {
	  //Ignore when we can't find the MWR devices
	  //   they don't always report and the timing of them can be annoying
	 
	}//catch
      }// Iterate over all the multiwire devices
    }// Iterate over all times
    

    //Here we will start collecting all the other beamline devices
    // First we get the times that the beamline device fired
    //  we have to pick a specific variable to use
    std::vector<double> times_temps = bfp->GetTimeList(fDeviceUsedForTiming);

    // We'll keep track of how many of these spills match to our 
    // DAQ trigger times
    int spill_count = 0;
    std::vector<int> matched_MWR;
    matched_MWR.resize(3);
    

    // Need to handle the first event in a run differently
    if(e.event() == 1){

      //We'll remove the spills after our event
      int spills_after_our_target = 0;
      // iterate through all the spills to find the 
      // spills that are after our triggered event
      for (size_t i = 0; i < times_temps.size(); i++) {       
	if(times_temps[i] > (t_current_event+fTimePad)){
	  spills_after_our_target++;
	}
      }//end loop through spill times 	 
      
      // Remove the spills after our trigger
      times_temps.erase(times_temps.end()-spills_after_our_target,times_temps.end());

      // Remove the spills before the start of our Run
      times_temps.erase(times_temps.begin(),times_temps.end()-number_of_gates_since_previous_event);

    }//end fix for "first event"

    // Iterating through each of the beamline times
    for (size_t i = 0; i < times_temps.size(); i++) {

      // Only continue if these times are matched to our DAQ time
      // plus or minus some time padding, currently using 3.3 ms 
      // which is half the Booster Rep Rate
      
      if(e.event() != 1){//We already addressed the "first event" above
	if(times_temps[i] > (t_current_event+fTimePad)){continue;}
	if(times_temps[i] <= (t_previous_event-fTimePad)){continue;}
      }

      //Loop through the multiwire devices:
      
      for(int dev = 0; dev < int(MWR_times.size()); dev++){

	//Loop through the multiwire times:
	double Tdiff = 1000000000.;
	matched_MWR[dev] = 0;
	bool best_match = false;
	for(int mwrt = 0;  mwrt < int(MWR_times[dev].size()); mwrt++){
	  
	  //found a candidate match! 
	  if(fabs((MWR_times[dev][mwrt] - times_temps[i])) < Tdiff){
	    best_match = true;
	    
	    //Check for a better match...
	    for (size_t j = 0; j < times_temps.size(); j++) {
	      if( j == i) continue;
	      if(times_temps[j] > (t_current_event+fTimePad)){continue;}
	      if(times_temps[j] <= (t_previous_event-fTimePad)){continue;}
	      
	      //is there a better match later in the spill sequence
	      if(fabs((MWR_times[dev][mwrt] - times_temps[j])) < 
		 fabs((MWR_times[dev][mwrt] - times_temps[i]))){
		//we can have patience...
		best_match = false;
		break;
	      }	     
	    }//end better match check

	    //Verified best match!
	    if(best_match == true){
	      matched_MWR[dev] = mwrt;
	      Tdiff = fabs((MWR_times[dev][mwrt] - times_temps[i]));
	    }
	  }//Find matches between 
	}//end loop over MWR times 
	
      }//end loop over MWR devices


      //Great we found a matched spill! Let's count it
      spill_count++;
      
      // initializing all of our device carriers
      // device definitions can be found in BNBSpillInfo.h

      double TOR860 = 0; // units e12 protons
      double TOR875 = 0; // units e12 protons
      double LM875A = 0; // units R/s
      double LM875B = 0; // units R/s
      double LM875C = 0; // units R/s
      double HP875 = 0; // units mm
      double VP875 = 0; // units mm
      double HPTG1 = 0; // units mm
      double VPTG1 = 0; // units mm
      double HPTG2 = 0; // units mm
      double VPTG2 = 0; // units mm
      double BTJT2 = 0; // units Deg C
      double THCURR = 0; // units kiloAmps
      
      double TOR860_time = 0; // units s
      
      // Here we request all the devices
      // since sometimes devices fail to report we'll
      // allow each to throw an exception but still move forward
      // interpreting these failures will be part of the beam quality analyses 
      try{bfp->GetNamedData(times_temps[i], "E:TOR860@",&TOR860,&TOR860_time);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:TOR875",&TOR875);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:LM875A",&LM875A);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:LM875B",&LM875B);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:LM875C",&LM875C);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:HP875",&HP875);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:VP875",&VP875);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:HPTG1",&HPTG1);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:VPTG1",&VPTG1);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:HPTG2",&HPTG2);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:VPTG2",&VPTG2);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:BTJT2",&BTJT2);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      try{bfp->GetNamedData(times_temps[i], "E:THCURR",&THCURR);}catch (WebAPIException &we) {std::cout << "At time : " << times_temps[i] << " " << "got exception: " << we.what() << "\n";}
      
      //crunch the times 
      unsigned long int time_closest_int = (int) TOR860_time;
      double time_closest_ns = (TOR860_time - time_closest_int)*1e9;
      
      //Store everything in our data-product
      sbn::BNBSpillInfo beamInfo;
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

      beamInfo.M875BB = unpacked_MWR[0][matched_MWR[0]];
      beamInfo.M875BB_spill_time_diff = (MWR_times[0][matched_MWR[0]] - times_temps[i]);
      beamInfo.M876BB = unpacked_MWR[1][matched_MWR[1]];
      beamInfo.M876BB_spill_time_diff = (MWR_times[1][matched_MWR[1]] - times_temps[i]);
      beamInfo.MMBTBB = unpacked_MWR[2][matched_MWR[2]];
      beamInfo.MMBTBB_spill_time_diff = (MWR_times[2][matched_MWR[2]] - times_temps[i]);

      fOutbeamInfos.push_back(beamInfo);
      // We do not write these to the art::Events because 
      // we can filter events but want to keep all the POT 
      // information, so we'll write it to the SubRun
    
    }//end iteration over beam device times

    if(spill_count > number_of_gates_since_previous_event)
      std::cout << "Event Spills : " << spill_count << ", DAQ Spills : " << number_of_gates_since_previous_event << " \t \t ::: WRONG!"<< std::endl;
    else
      std::cout << "Event Spills : " << spill_count << ", DAQ Spills : " << number_of_gates_since_previous_event << std::endl;

  } //end check if BNB DAQ triggered gate
}//end iteration over art::Events

void sbn::BNBRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

//____________________________________________________________________________                                                                                                                                                                                      
void sbn::BNBRetriever::endSubRun(art::SubRun& sr)
{
  // We will add all of the BNBSpillInfo data-products to the 
  // art::SubRun so it persists 
  // currently this is ~2.7 kB/event or ~0.07 kB/spill

std::cout << "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
std::cout << "Total number of Selected Spills : " << fOutbeamInfos.size() << std::endl;

  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >(fOutbeamInfos);

  sr.put(std::move(p));

  return;
}

DEFINE_ART_MODULE(sbn::BNBRetriever)    
