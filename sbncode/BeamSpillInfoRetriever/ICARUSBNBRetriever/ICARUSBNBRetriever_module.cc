/** ********************************************************************
 * @file ICARUSBNBRetriever_module.cc
 * @date Wed April 9 2021
 * @author J. Zennamo (FNAL)
 * 
 * Based heavily on code by Z. Pavlovic written for MicroBooNE 
 * Based heavily on code by NOvA collaboration (Thanks NOvA!):
 *        https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/BNBInfo_module.cc
 * Database implementation by Justin Mueller 
 */

#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbncode/BeamSpillInfoRetriever/SBNDPOTTools.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbndaq-artdaq-core/Overlays/ICARUS/ICARUSTriggerV3Fragment.hh"
#include <sqlite3.h>

namespace sbn {
  class ICARUSBNBRetriever;
}

class sbn::ICARUSBNBRetriever : public art::EDProducer {
public:
  
  struct Config {
    
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    
    fhicl::Atom<double> TimePadding {
      Name{ "TimePadding" },
      Comment{ "extension to the time window considered when collecting spills [seconds]" },
      0.0333 // default
      };
 
    fhicl::Atom<double> BESOffset {
      Name{ "BESOffset" },
      Comment{ "Offset between beam early warning signal and $1D [seconds]" },
      0.036 // default
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

    fhicl::Atom<std::string> TriggerDatabaseFile {
      Name{ "TriggerDatabaseFile" },
      Comment{ "" } // explain what this is, what's for and its unit
      };
    

  }; // Config
  
  using Parameters = art::EDProducer::Table<Config>;
  
  
  explicit ICARUSBNBRetriever(Parameters const& params);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ICARUSBNBRetriever(ICARUSBNBRetriever const&) = delete;
  ICARUSBNBRetriever(ICARUSBNBRetriever&&) = delete;
  ICARUSBNBRetriever& operator=(ICARUSBNBRetriever const&) = delete;
  ICARUSBNBRetriever& operator=(ICARUSBNBRetriever&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

private:
  // input labels
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  double fTimePad;
  double fBESOffset;
  std::string fURL;
  MWRData mwrdata;
  int run_number;
  std::string raw_data_label;
  std::string fDeviceUsedForTiming;
  unsigned int TotalBeamSpills;  
  //
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;
  
  //
  std::string fTriggerDatabaseFile;
  sqlite3 *db;
  int rc;

  static constexpr double MWRtoroidDelay = -0.035; ///< the same time point is measured _t_ by MWR and _t + MWRtoroidDelay`_ by the toroid [ms]

  /// Returns the information of the trigger in the current event.
  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  /**
   * @brief Matches spill times with multiwire chamber data from the database.
   * @param eventID ID of the event the information is associated to
   * @param triggerInfo information from the trigger of this event
   * @param MWRdata data from multiwire chambers
   * @param isFirstEventInRun whether we are processing the first event of the run
   * @param[out] beamInfos container to _add_ spill information records to
   * @return count of matched spills
   */
  int matchMultiWireData(
    art::EventID const& eventID, 
    TriggerInfo_t const& triggerInfo,
    MWRdata_t const& MWRdata, bool isFirstEventInRun,
    std::vector< sbn::BNBSpillInfo >& beamInfos
    ) const;

/**
 * @brief SQLite callback function for retrieving trigger_type from a query.
 * @param data Pointer to the integer where the trigger type will be stored.
 * @param argc Count of the number of columns returned by the query.
 * @param argv Array of c-strings containing the column data of the query.
 * @param columns Array of c-strings listing the names of the columns.
 * @return 0 if successful.
 */
//  int callback_trigger_type(void *data, 	
//			    int argc, 
//			    char **argv, 
//			    char **columns); 
  
/**
 * @brief Queries the trigger database and finds the trigger_type of the matching trigger (if any).
 * @param db The pointer to the SQLite database instance.
 * @param run The run number of the current event (helps with queries).
 * @param gate_time The time in milliseconds of the gate.
 * @param threshold The required absolute time difference between gate and trigger.
 * @return trigger_type -1: No matching trigger, 0: Majority, 1: MinBias
 */
  int get_trigger_type_matching_gate(sqlite3 *db, 
				     int func(void*,int,char**,char**), 
				     int run, 
				     long long int gate_time, 
				     float threshold) const;
  
};

int callback_trigger_type(void *data, int argc, char **argv, char **columns)
{
  int *result = static_cast<int*>(data);
  // Does this query return non-NULL values?
  if(argc > 0 && argv[0])
    *result = std::stoi(argv[0]);
  else
    *result = -1;

  return 0;
}

int sbn::ICARUSBNBRetriever::get_trigger_type_matching_gate(sqlite3 *db, int func(void*,int,char**,char**), int run, long long int gate_time, float threshold) const
{
  int trigger_type(-1), query_status;
  std::stringstream query;
  query << "SELECT trigger_type FROM triggerdata WHERE gate_type=1 AND run_number ="
        << run
        << " AND ABS(1000000000*wr_seconds + wr_nanoseconds  - "
        << std::fixed << gate_time
        << ") < "
        << threshold*1000000
        << " ORDER BY ABS(1000000000*wr_seconds + wr_nanoseconds  - "
	<< std::fixed << gate_time
	<< ") LIMIT 1;";

  query_status = sqlite3_exec(db, query.str().c_str(), func, &trigger_type, NULL);
  if (query_status != SQLITE_OK)
  {
    mf::LogError("BNBEXTRetriever") << "SQL error: " << sqlite3_errmsg(db);
    trigger_type = -1;
  }
  return trigger_type;
}

sbn::ICARUSBNBRetriever::ICARUSBNBRetriever(Parameters const& params)
  : EDProducer{params},
  fTimePad(params().TimePadding()),
  fBESOffset(params().BESOffset()),
  raw_data_label(params().RawDataLabel()),
  fDeviceUsedForTiming(params().DeviceUsedForTiming()),
  bfp(     ifbeam_handle->getBeamFolder(params().Bundle(), params().URL(), params().TimeWindow())),
  bfp_mwr( ifbeam_handle->getBeamFolder(params().MultiWireBundle(), params().URL(), params().MWR_TimeWindow())),
  fTriggerDatabaseFile(params().TriggerDatabaseFile())
{
  
  // Check fTimePad is positive 
  if (fTimePad < 0) {
    throw art::Exception(art::errors::Configuration)
      << "Parameter `TimePadding` must be non-negative (" << fTimePad << " was specified).\n";
  }//End Time check  
  
  // how close in time does the spill time have to be from the DAQ time (in seconds).
  // If these are too large then it fails to capture the device 
  // If these are too small then the time jitter in devices means we miss good data 
  //
  // These values should likely not be changed unless authors of the IFBeam API are consulted
  //
  bfp->set_epsilon(0.02); //20 ms, this was tuned by hand and compared to IFBeamDB times  
  bfp_mwr->set_epsilon(0.5);

  //bfp_mwr->setValidWindow(86400);  
  bfp_mwr->setValidWindow(3605);  
  produces< std::vector< sbn::BNBSpillInfo >, art::InSubRun >();
  TotalBeamSpills = 0;

  cet::search_path sp("FW_SEARCH_PATH");
  std::string trigDB_path = sp.find_file(fTriggerDatabaseFile.c_str());

  rc = sqlite3_open(trigDB_path.c_str(), &db);
  if(rc)
    {
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      throw art::Exception(art::errors::NotFound)
	<< "Can't open database: " << sqlite3_errmsg(db);
    }

}


void sbn::ICARUSBNBRetriever::produce(art::Event& e)
{

  // If this is the first event in the run, then ignore it
  // We do not currently have the ability to figure out the first
  // spill that the DAQ was sensitive to, so don't try to save any
  // spill information
  //
  // TODO: long-term goal -- can we fix this?
  // FIXME This is wrong.... 
  //      Need to use:  ICARUSTriggerV3Fragment long getTotalTriggerBNBMaj() const

  if (e.event() == 1) return;
  
  run_number = e.id().run();

  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);
  
  //We only want to process BNB gates, i.e. type 1 
  if(triggerInfo.gate_type != 1) return;
  // Keep track of the number of beam gates the DAQ thinks 
  //   are in this job
  TotalBeamSpills += triggerInfo.number_of_gates_since_previous_event;
  MWRdata_t const MWRdata = extractSpillTimes(triggerInfo, bfp, bfp_mwr, fTimePad, MWRtoroidDelay, mwrdata);
  
  int const spill_count = matchMultiWireData(e.id(), triggerInfo, MWRdata, e.event() == 1, fOutbeamInfos);
  
  if(spill_count > int(triggerInfo.number_of_gates_since_previous_event))
    mf::LogDebug("ICARUSBNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << " \t \t ::: WRONG!"<< std::endl;
  else
    mf::LogDebug("ICARUSBNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << std::endl;
  
}//end iteration over art::Events


sbn::TriggerInfo_t sbn::ICARUSBNBRetriever::extractTriggerInfo(art::Event const& e) const {
  
  //Here we read in the artdaq Fragments and extract three pieces of information:
  // 1. The time of the current event, t_current_event
  // 2. the time of the previously triggered event, t_previous_event (NOTE: Events are non-sequential!)
  // 3. the number of beam spills since the previously triggered event, number_of_gates_since_previous_event
  
  auto const & raw_data = e.getProduct< std::vector<artdaq::Fragment> >({ raw_data_label, "ICARUSTriggerV3" });
  auto const & extraTrigInfo = e.getProduct< sbn::ExtraTriggerInfo >("daqTrigger");

  TriggerInfo_t triggerInfo;
  
  triggerInfo.WR_to_Spill_conversion = extraTrigInfo.WRtimeToTriggerTime;   
  
  for(auto raw_datum : raw_data){
   
    uint64_t artdaq_ts = raw_datum.timestamp();
    icarus::ICARUSTriggerV3Fragment frag(raw_datum);
    std::string data = frag.GetDataString();
    char *buffer = const_cast<char*>(data.c_str());
    icarus::ICARUSTriggerInfo datastream_info = icarus::parse_ICARUSTriggerV3String(buffer);
    triggerInfo.gate_type = datastream_info.gate_type;
    triggerInfo.number_of_gates_since_previous_event = frag.getDeltaGatesBNBMaj();
    
    /*                                                                                                                  
       The DAQ trigger time is issued at the Beam Extraction Signal (BES) which is issued                               
       36 ms *after* the $1D of the BNB, which is what is used in the IFBeam database                                   
                                                                                                                        
       We subtract 36ms from the Trigger time to match our triggers to the spills in the                                
       IFBeam database                                                                                                  
                                                                                                                        
    */

    double BESinTSUnits = fBESOffset*1e9; // flc param is in seconds need to convert to match TS
    triggerInfo.t_current_event = static_cast<double>(artdaq_ts-BESinTSUnits)/(1000000000.0); //check this offset...
    if(triggerInfo.gate_type == 1)
      triggerInfo.t_previous_event = (static_cast<double>(frag.getLastTimestampBNBMaj()-BESinTSUnits))/(1e9);
    else
      triggerInfo.t_previous_event = (static_cast<double>(frag.getLastTimestampOther()-BESinTSUnits))/(1000000000.0);
    
  }
  
  mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "Previous : " << triggerInfo.t_previous_event << ", Current : " << triggerInfo.t_current_event << ", Spill Count " << triggerInfo.number_of_gates_since_previous_event << std::endl;

  return triggerInfo;
}

int sbn::ICARUSBNBRetriever::matchMultiWireData(
  art::EventID const& eventID,
  TriggerInfo_t const& triggerInfo,
  MWRdata_t const& MWRdata, bool isFirstEventInRun,
  std::vector< sbn::BNBSpillInfo >& beamInfos
) const {
  
  auto const& [ MWR_times, unpacked_MWR ] = MWRdata; // alias
  
  //Here we will start collecting all the other beamline devices
  // First we get the times that the beamline device fired
  //  we have to pick a specific variable to use
  std::vector<double> times_temps = bfp->GetTimeList(fDeviceUsedForTiming);
  
  mf::LogDebug("ICARUSBNBRetriever") << "matchMultiWireData:: Number of time spills : " << times_temps.size() << std::endl;

  // We'll keep track of how many of these spills match to our 
  // DAQ trigger times
  int spill_count = 0;
  int spills_removed = 0;
  std::vector<int> matched_MWR;
  matched_MWR.resize(3);
  
  
  // NOTE: for now, this is dead code because we don't
  // do anything for the first event in a run. We may want to revisit 
  // this later to understand if there is a way we can do the POT
  // accounting in the first event.
  //
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
    times_temps.erase(times_temps.begin(), times_temps.end() - std::min(int(triggerInfo.number_of_gates_since_previous_event), int(times_temps.size())));
        
  }//end fix for "first event"
  
    ///reject time_stamps which have a trigger_type == 1 from data-base
    //To-Do 

  //  mf::LogDebug("ICARUSBNBRetriever") << "Total number of Times we're going to test: " << times_temps.size() <<  std::endl;
  // mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "Upper Limit : " << (triggerInfo.t_current_event)+fTimePad <<  std::endl;
  // mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "Lower Limit : " << (triggerInfo.t_previous_event)+fTimePad <<  std::endl;
  
  // Iterating through each of the beamline times
  for (size_t i = 0; i < times_temps.size(); i++) {
    
    // Only continue if these times are matched to our DAQ time
    //mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "Time # : " <<  i << std::endl;

    if(!isFirstEventInRun){//We already addressed the "first event" above
      if(times_temps[i] > (triggerInfo.t_current_event)+fTimePad){
	//mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "Removed!  : " << times_temps[i] << std::endl;
	spills_removed++; 
	continue;} 
      if(times_temps[i] <= (triggerInfo.t_previous_event)+fTimePad){
	spills_removed++; 
	//mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "Removed!  : " << times_temps[i] << std::endl;
	continue;}
    }

    //check if this spill is is minbias   
    /*
      40 ms was selected to be close to but outside the 66 ms 
      time of the next spill (when the beam is running at 15 Hz) 
      DocDB 33155 provides documentation of this
    */

    double BESinTSUnits = fBESOffset*1e9; // flc param is in seconds need to convert to match TS
    mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19) << "matchMultiWireData:: trigger type : " << get_trigger_type_matching_gate(db, callback_trigger_type, run_number, times_temps[i]*1.e9-triggerInfo.WR_to_Spill_conversion+BESinTSUnits, 40.) << " times : spill " << times_temps[i]*1.e9 << " - " << triggerInfo.WR_to_Spill_conversion << " + " << BESinTSUnits <<  std::endl;
    
    if(get_trigger_type_matching_gate(db, callback_trigger_type, run_number, times_temps[i]*1.e9-triggerInfo.WR_to_Spill_conversion+BESinTSUnits, 40.) == 1){
          mf::LogDebug("ICARUSBNBRetriever") << std::setprecision(19)  << "matchMultiWireData:: Skipped a MinBias gate at : " << times_temps[i]*1000. << std::endl;

      continue;
    }
      
    //Great we found a matched spill! Let's count it
    spill_count++;

    //Loop through the multiwire devices:
    
    for(int dev = 0; dev < int(MWR_times.size()); dev++){
      
      //Loop through the multiwire times:
      double Tdiff = 1000000000.;
      matched_MWR[dev] = 0;

      for(int mwrt = 0;  mwrt < int(MWR_times[dev].size()); mwrt++){

	//found a candidate match! 
	if(fabs((MWR_times[dev][mwrt] - times_temps[i])) >= Tdiff){continue;}
	
	bool best_match = true;
	  
	//Check for a better match...
	for (size_t j = 0; j < times_temps.size(); j++) {
	  if( j == i) continue;
	  if(times_temps[j] > (triggerInfo.t_current_event+fTimePad)){continue;}
	  if(times_temps[j] <= (triggerInfo.t_previous_event+fTimePad)){continue;}
	  
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
	
      }//end loop over MWR times 
      
    }//end loop over MWR devices
    
    sbn::BNBSpillInfo spillInfo = makeBNBSpillInfo(eventID, times_temps[i], MWRdata, matched_MWR, bfp);

    beamInfos.push_back(std::move(spillInfo));

    // We do not write these to the art::Events because 
    // we can filter events but want to keep all the POT 
    // information, so we'll write it to the SubRun
    
  }//end iteration over beam device times
  
  //  mf::LogDebug("ICARUSBNBRetriever") << "matchMultiWireData:: Total spills counted:  " << spill_count << "   Total spills removed : " << spills_removed <<  std::endl;

  return spill_count;
}

void sbn::ICARUSBNBRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

//____________________________________________________________________________                                                                                                                                                                                      
void sbn::ICARUSBNBRetriever::endSubRun(art::SubRun& sr)
{
  // We will add all of the BNBSpillInfo data-products to the 
  // art::SubRun so it persists 
  // currently this is ~2.7 kB/event or ~0.07 kB/spill

mf::LogDebug("ICARUSBNBRetriever")<< "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
mf::LogDebug("ICARUSBNBRetriever")<< "Total number of Selected Spills : " << fOutbeamInfos.size() << std::endl;

  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  std::swap(*p, fOutbeamInfos);

  sr.put(std::move(p), art::subRunFragment());

  return;
}

DEFINE_ART_MODULE(sbn::ICARUSBNBRetriever)    
