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
#include "fhiclcpp/types/Atom.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/Utilities/sparse_vector.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcorealg/CoreUtils/counter.h"
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
#include <vector>

namespace sbn {
  class BNBRetriever;
}

class sbn::BNBRetriever : public art::EDProducer {
public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<double> TimePadding {
      Name{ "TimePadding" },
      Comment{ "extension to the time window considered when collecting spills [seconds]" },
      0.0333 // default
      };
    
    fhicl::Atom<std::string> RawDataLabel {
      Name{ "raw_data_label" },
      Comment{ "art data product instance name for trigger information (product label is 'daq')" }
      };
    
    fhicl::Atom<std::string> DeviceUsedForTiming {
      Name{ "DeviceUsedForTiming" },
      Comment{ "name in the IFBeam database of the device used to extract spill times" }
      };
    
    fhicl::Atom<std::string> URL {
      Name{ "URL" },
      Comment{ "IFBeam database access URL" }
      };
    
    fhicl::Atom<std::string> Bundle {
      Name{ "Bundle" },
      Comment{ "" } // explain what this is and which database/table it's looking for
      };
    
    fhicl::Atom<double> TimeWindow {
      Name{ "TimeWindow" },
      Comment{ "" } // explain what this is, what's for and its unit
      };
    
    fhicl::Atom<std::string> MultiWireBundle {
      Name{ "MultiWireBundle" },
      Comment{ "" } // explain what this is and which database/table it's looking for
      };
    
    fhicl::Atom<double> MWR_TimeWindow {
      Name{ "MWR_TimeWindow" },
      Comment{ "" } // explain what this is, what's for and its unit
      };
    
  }; // Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  explicit BNBRetriever(Parameters const& params);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BNBRetriever(BNBRetriever const&) = delete;
  BNBRetriever(BNBRetriever&&) = delete;
  BNBRetriever& operator=(BNBRetriever const&) = delete;
  BNBRetriever& operator=(BNBRetriever&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

private:
  // input labels
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  double fTimePad;
  std::string fURL;
  MWRData mwrdata;
  std::string raw_data_label_;
  std::string fDeviceUsedForTiming;
  unsigned int TotalBeamSpills;  
  //
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;

  struct TriggerInfo_t {
    int gate_type = 0; ///< Source of the spill: `1`: BNB, `2`: NuMI
    double t_current_event  = 0;
    double t_previous_event = 0;
    double number_of_gates_since_previous_event = 0; // FIXME needs to be integral type
  };
  
  struct MWRdata_t {
    std::vector< std::vector<double> > MWR_times;
    std::vector< std::vector< std::vector< int > > > unpacked_MWR;
  };
  

  static constexpr double MWRtoroidDelay = -0.035; ///< the same time point is measured _t_ by MWR and _t + MWRtoroidDelay`_ by the toroid [ms]

  /// Returns the information of the trigger in the current event.
  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  
  /**
   * @brief Determines spill times and extracts data based on multiwire devices.
   * @param triggerInfo information from the trigger of this event
   * @return times and unpacked data, per device (`"M875BB"`, `"M876BB"`, `"MMBTBB"`)
   */
  MWRdata_t extractSpillTimes(TriggerInfo_t const& triggerInfo) const;

  /**
   * @brief Matches spill times with multiwire chamber data from the database.
   * @param triggerInfo information from the trigger of this event
   * @param MWRdata data from multiwire chambers
   * @param isFirstEventInRun whether we are processing the first event of the run
   * @param[out] beamInfos container to _add_ spill information records to
   * @return count of matched spills
   */
  int matchMultiWireData(
    TriggerInfo_t const& triggerInfo,
    MWRdata_t const& MWRdata, bool isFirstEventInRun,
    std::vector< sbn::BNBSpillInfo >& beamInfos
    ) const;

  /**
   * @brief Assembles and returns a spill information record.
   * @param time time of the spill
   * @param MWRdata all extracted data from multiwire chambers
   * @param matched_MWR data from multiwire chambers matched with the time
   * @return a `sbn::BNBSpillInfo` object with information on the spill at `time`
   */
  sbn::BNBSpillInfo makeBNBSpillInfo
    (double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR) const;

};

sbn::BNBRetriever::BNBRetriever(Parameters const& params)
  : EDProducer{params},
  fTimePad(params().TimePadding()),
  raw_data_label_(params().RawDataLabel()),
  fDeviceUsedForTiming(params().DeviceUsedForTiming()),
  bfp(     ifbeam_handle->getBeamFolder(params().Bundle(), params().URL(), params().TimeWindow())),
  bfp_mwr( ifbeam_handle->getBeamFolder(params().MultiWireBundle(), params().URL(), params().MWR_TimeWindow()))
{
  
  // Check fTimePad is positive 
  if (fTimePad < 0) {
    throw art::Exception(art::errors::Configuration)
      << "Parameter `TimePadding` must be non-negative (" << fTimePad << " was specified).\n";
  }//End Time check  
  
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
  
  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);
  
  //We only want to process BNB gates, i.e. type 1 
  if(triggerInfo.gate_type != 1) return;
  // Keep track of the number of beam gates the DAQ thinks 
  //   are in this job
  TotalBeamSpills += triggerInfo.number_of_gates_since_previous_event;
  
  
  MWRdata_t const MWRdata = extractSpillTimes(triggerInfo);
  
  
  int const spill_count = matchMultiWireData(triggerInfo, MWRdata, e.event() == 1, fOutbeamInfos);
  
  
  if(spill_count > triggerInfo.number_of_gates_since_previous_event)
    mf::LogDebug("BNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << " \t \t ::: WRONG!"<< std::endl;
    else
      mf::LogDebug("BNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << std::endl;
  
}//end iteration over art::Events


sbn::BNBRetriever::TriggerInfo_t sbn::BNBRetriever::extractTriggerInfo(art::Event const& e) const {
  
  //Here we read in the artdaq Fragments and extract three pieces of information:
  // 1. The time of the current event, t_current_event
  // 2. the time of the previously triggered event, t_previous_event (NOTE: Events are non-sequential!)
  // 3. the number of beam spills since the previously triggered event, number_of_gates_since_previous_event
  
  auto const & raw_data = e.getByLabel< std::vector<artdaq::Fragment> >({ raw_data_label_, "ICARUSTriggerUDP" });
  
  TriggerInfo_t triggerInfo;

  for(auto raw_datum : raw_data){
   
    uint64_t artdaq_ts = raw_datum.timestamp();
    icarus::ICARUSTriggerUDPFragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerString(buffer);
    triggerInfo.gate_type = datastream_info.gate_type;
    triggerInfo.number_of_gates_since_previous_event = frag.getDeltaGatesBNB();
  
    triggerInfo.t_current_event = static_cast<double>(artdaq_ts)/(1000000000.0); //check this offset...
    if(triggerInfo.gate_type == 1)
      triggerInfo.t_previous_event = (static_cast<double>(frag.getLastTimestampBNB()))/(1e9);
    else
      triggerInfo.t_previous_event = (static_cast<double>(frag.getLastTimestampOther()))/(1000000000.0);
    
  }
  
  mf::LogDebug("BNBRetriever") << std::setprecision(19) << "Previous : " << triggerInfo.t_previous_event << ", Current : " << triggerInfo.t_current_event << std::endl;

  return triggerInfo;
}


sbn::BNBRetriever::MWRdata_t sbn::BNBRetriever::extractSpillTimes(TriggerInfo_t const& triggerInfo) const {
  
  // These lines get everything primed within the IFBeamDB
  //   They seem redundant but they are needed
  try{auto cur_vec_temp = bfp->GetNamedVector((triggerInfo.t_previous_event)-fTimePad,"E:THCURR");} catch (WebAPIException &we) {}      
  try{auto packed_M876BB_temp = bfp_mwr->GetNamedVector((triggerInfo.t_current_event)+fTimePad,"E:M875BB{4440:888}.RAW");} catch (WebAPIException &we) {}
  
  
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

  int t_steps = int(((triggerInfo.t_previous_event - fTimePad) - (triggerInfo.t_current_event + fTimePad))/0.5)+25;
  
  for(int t = 0; t < t_steps; t++){//Iterate through time increments
    for (std::string const& var : vars) {// Iterate through the devices
      
      //Make sure we have a device
      if(var.empty()) continue;
      
      /// Check the device name and interate the double-vector index
      if(var.find("M875BB") != std::string::npos ) dev = 0;
      else if(var.find("M876BB") != std::string::npos ) dev = 1;
      else if(var.find("MMBTBB") != std::string::npos ) dev = 2;
      else{continue;}
      
      time_for_mwr = 0;
      
      try{
	//Pull the MWR data for the device
	// these data are "packed"
	std::vector<double> packed_MWR = bfp_mwr->GetNamedVector((triggerInfo.t_previous_event)-fTimePad+double(0.5*t),var,&time_for_mwr);

	//We'll convert this into a format
	// that we can unpack doubles >> strings
	//
	packed_data_str.clear();
	packed_data_str += std::to_string(int(time_for_mwr));
	packed_data_str.append(",");
	packed_data_str.append(var);
	packed_data_str.append(",,");
	
	/*	for(auto const value: packed_MWR){
	  packed_data_str += ',';
	  packed_data_str += std::to_string(int(value));
	  }*/
	for(int j = 0; j < int(packed_MWR.size()); j++){
	  packed_data_str += std::to_string(int(packed_MWR[j]));
	  if(j < int(packed_MWR.size())-1)
	    packed_data_str.append(",");
	}
	
	// Use Zarko's unpacking function to turn this into consumeable data
	std::vector<double> MWR_times_temp;
	
	// There is a 35 ms offset between the toriod and the MWR times
	//   we'll just remove that here to match to the spill times
	std::vector< std::vector< int > > unpacked_MWR_temp = mwrdata.unpackMWR(packed_data_str,MWR_times_temp,MWRtoroidDelay);
	
	//There are four events that are packed into one MWR IFBeam entry
	for(std::size_t s: util::counter(unpacked_MWR_temp.size())){
	  	  
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
  
  return { std::move(MWR_times), std::move(unpacked_MWR) };
}


int sbn::BNBRetriever::matchMultiWireData(
  TriggerInfo_t const& triggerInfo,
  MWRdata_t const& MWRdata, bool isFirstEventInRun,
  std::vector< sbn::BNBSpillInfo >& beamInfos
) const {
  
  auto const& [ MWR_times, unpacked_MWR ] = MWRdata; // alias
  
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
  if(isFirstEventInRun){
    
    //We'll remove the spills after our event
    int spills_after_our_target = 0;
    // iterate through all the spills to find the 
    // spills that are after our triggered event
    for (size_t i = 0; i < times_temps.size(); i++) {       
      if(times_temps[i] > (triggerInfo.t_current_event+fTimePad)){
	spills_after_our_target++;
      }
    }//end loop through spill times 	 
    
    // Remove the spills after our trigger
    times_temps.erase(times_temps.end()-spills_after_our_target,times_temps.end());
    
    // Remove the spills before the start of our Run
    times_temps.erase(times_temps.begin(),times_temps.end()-triggerInfo.number_of_gates_since_previous_event);
    
  }//end fix for "first event"
  
  // Iterating through each of the beamline times
  for (size_t i = 0; i < times_temps.size(); i++) {
    
    // Only continue if these times are matched to our DAQ time
    // plus or minus some time padding, currently using 3.3 ms 
    // which is half the Booster Rep Rate
    
    if(!isFirstEventInRun){//We already addressed the "first event" above
      if(times_temps[i] > (triggerInfo.t_current_event+fTimePad)){continue;}
      if(times_temps[i] <= (triggerInfo.t_previous_event-fTimePad)){continue;}
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
	    if(times_temps[j] > (triggerInfo.t_current_event+fTimePad)){continue;}
	    if(times_temps[j] <= (triggerInfo.t_previous_event-fTimePad)){continue;}
	    
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
    
    sbn::BNBSpillInfo spillInfo = makeBNBSpillInfo(times_temps[i], MWRdata, matched_MWR);

    beamInfos.push_back(std::move(spillInfo));
    // We do not write these to the art::Events because 
    // we can filter events but want to keep all the POT 
    // information, so we'll write it to the SubRun
    
  }//end iteration over beam device times
  
  return spill_count;
}


sbn::BNBSpillInfo sbn::BNBRetriever::makeBNBSpillInfo
  (double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR) const
{
  
  auto const& [ MWR_times, unpacked_MWR ] = MWRdata; // alias
  
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
  try{bfp->GetNamedData(time, "E:TOR860@",&TOR860,&TOR860_time);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:TOR875",&TOR875);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:LM875A",&LM875A);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:LM875B",&LM875B);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:LM875C",&LM875C);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:HP875",&HP875);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:VP875",&VP875);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:HPTG1",&HPTG1);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:VPTG1",&VPTG1);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:HPTG2",&HPTG2);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:VPTG2",&VPTG2);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:BTJT2",&BTJT2);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:THCURR",&THCURR);}catch (WebAPIException &we) {mf::LogDebug("BNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  
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
  beamInfo.M875BB_spill_time_diff = (MWR_times[0][matched_MWR[0]] - time);
  beamInfo.M876BB = unpacked_MWR[1][matched_MWR[1]];
  beamInfo.M876BB_spill_time_diff = (MWR_times[1][matched_MWR[1]] - time);
  beamInfo.MMBTBB = unpacked_MWR[2][matched_MWR[2]];
  beamInfo.MMBTBB_spill_time_diff = (MWR_times[2][matched_MWR[2]] - time);
  
  // We do not write these to the art::Events because 
  // we can filter events but want to keep all the POT 
  // information, so we'll write it to the SubRun
  
  return beamInfo;
}


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

mf::LogDebug("BNBRetriever")<< "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
mf::LogDebug("BNBRetriever")<< "Total number of Selected Spills : " << fOutbeamInfos.size() << std::endl;

  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >(fOutbeamInfos);

  sr.put(std::move(p), art::subRunFragment());

  return;
}

DEFINE_ART_MODULE(sbn::BNBRetriever)    
