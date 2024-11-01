////////////////////////////////////////////////////////////////////////
// Class:       SBNDBNBRetriever
// Plugin Type: producer 
// File:        SBNDBNBRetriever_module.cc
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

#include <memory>
#include <bitset>
#include <tuple>
#include <algorithm>

#include "sbndaq-artdaq-core/Overlays/SBND/PTBFragment.hh"
#include "sbndaq-artdaq-core/Overlays/SBND/TDCTimestampFragment.hh"
#include "artdaq-core/Data/ContainerFragment.hh"
//#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"

#include "ifdh_art/IFBeamService/IFBeam_service.h"
#include "ifbeam_c.h"
#include "MWRData.h"

#include "larcorealg/CoreUtils/counter.h"

class SBNDBNBRetriever : public art::EDProducer {
public:
  explicit SBNDBNBRetriever(fhicl::ParameterSet const & params);
  
  // Plugins should not be copied or assigned.
  SBNDBNBRetriever(SBNDBNBRetriever const &) = delete;
  SBNDBNBRetriever(SBNDBNBRetriever &&) = delete;
  SBNDBNBRetriever & operator = (SBNDBNBRetriever const &) = delete;
  SBNDBNBRetriever & operator = (SBNDBNBRetriever &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

private:
  // Declare member data here.
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  double fTimePad;
  std::string fInputLabel;
  std::string fInputNonContainerInstance;
  std::string fDeviceUsedForTiming;
  std::string fOutputInstance;
  std::string raw_data_label;
  int fDebugLevel;
  sbn::MWRData mwrdata;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;

  struct TriggerInfo_t {
    double t_current_event  = 0;
    double t_previous_event = 0;
    unsigned int number_of_gates_since_previous_event = 0; // FIXME needs to be integral type
  };

  struct MWRdata_t {
    std::vector< std::vector<double> > MWR_times;
    std::vector< std::vector< std::vector< int > > > unpacked_MWR;
  };

  static constexpr double MWRtoroidDelay = -0.035; ///< the same time point is measured _t_ by MWR and _t + MWRtoroidDelay`_ by the toroid [ms]

  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  double extractPTBTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const;
  double extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const;
  MWRdata_t extractSpillTimes(TriggerInfo_t const& triggerInfo) const; 
  int matchMultiWireData(
    art::EventID const& eventID, 
    TriggerInfo_t const& triggerInfo,
    MWRdata_t const& MWRdata, bool isFirstEventInRun,
    std::vector< sbn::BNBSpillInfo >& beamInfos
    ) const;

  sbn::BNBSpillInfo makeBNBSpillInfo
    (art::EventID const& eventID, double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR) const;
};

SBNDBNBRetriever::SBNDBNBRetriever(fhicl::ParameterSet const & params)
  : EDProducer{params} {
  raw_data_label = params.get<std::string>("raw_data_label", "daq");
  fInputLabel = params.get<std::string>("InputLabel");
  fDeviceUsedForTiming = params.get<std::string>("DeviceUsedForTiming");
  fTimePad = params.get<double>("TimePadding");
  fInputNonContainerInstance = params.get<std::string>("InputNonContainerInstance");
  fOutputInstance = params.get<std::string>("OutputInstance");
  fDebugLevel = params.get<int>("DebugLevel",0);
  bfp = ifbeam_handle->getBeamFolder(params.get<std::string>("Bundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("TimeWindow")));
  //bfp->set_epsilon(0.02);
  bfp->set_epsilon(0.05);
  bfp_mwr = ifbeam_handle->getBeamFolder(params.get<std::string>("MultiWireBundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("MWR_TimeWindow")));
  bfp_mwr->set_epsilon(0.5);
  bfp_mwr->setValidWindow(3605);
}

int eventNum =0;
int _run;
int _subrun;
int _event;

void SBNDBNBRetriever::produce(art::Event & e)
{

  // If this is the first event in the run, then ignore it
  // We do not currently have the ability to figure out the first
  // spill that the DAQ was sensitive to, so don't try to save any
  // spill information

  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);

  MWRdata_t const MWRdata = extractSpillTimes(triggerInfo);

  //std::cout << "event: " << e.event() << std::endl;
  int const spill_count = matchMultiWireData(e.id(), triggerInfo, MWRdata, e.event() == 1, fOutbeamInfos);

  if(spill_count > int(triggerInfo.number_of_gates_since_previous_event))
    mf::LogDebug("SBNDBNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << " \t \t ::: WRONG!"<< std::endl;
  else
    mf::LogDebug("SBNDBNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << std::endl;
}

double SBNDBNBRetriever::extractPTBTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const {
  int numcont = 0;
  uint64_t PTBTimeStamp = 0;
  uint64_t PTBTriggerWord;
  for (auto const& cont : *cont_frags)
  { 
    artdaq::ContainerFragment cont_frag(cont);
    numcont++;
    int numfrag = 0;
    for (size_t fragi = 0; fragi < cont_frag.block_count(); ++fragi)
    {
      numfrag++;
      artdaq::Fragment frag = *cont_frag[fragi];
      sbndaq::CTBFragment ctb_frag(frag);   // somehow the name CTBFragment stuck
      for(size_t word_i = 0; word_i < ctb_frag.NWords(); ++word_i)
      {
        if(ctb_frag.Trigger(word_i)){
          uint32_t wt = 0;
          uint32_t word_type = ctb_frag.Word(word_i)->word_type;
          wt = word_type;
	  PTBTriggerWord = ctb_frag.Trigger(word_i)->trigger_word & 0x1FFFFFFFFFFFFFFF; //& 0x1FFFFFFFFFFF extracts the 61 bit payload
          std::bitset<32> desired_word{"00000000000000000000000000000010"};// Only use timestamp from specific HLT trigger word  
	  if (wt == 2 && std::bitset<32>(PTBTriggerWord) == desired_word)
	  {
	    PTBTimeStamp = ctb_frag.Trigger(word_i)->timestamp * 20; 
	    std::cout << "PTB Word type [" << std::bitset<3>(ctb_frag.Word(word_i)->word_type) << "]  ";
	    std::cout << std::bitset<32>(PTBTriggerWord) << "  HLT Timestamp: "  << std::bitset<64>(PTBTimeStamp/20) <<std::endl;     
	  }
        }
      } //End of loop over the number of trigger words
    } //End of loop over the number of fragments per container
  } //End of loop over the number of containers
  return std::bitset<64>(PTBTimeStamp/20).to_ullong()/50e6; 
}

double SBNDBNBRetriever::extractTDCTimeStamp(art::Handle<std::vector<artdaq::Fragment> > cont_frags) const {
  int numcont = 0;
  uint64_t TDCTimeStamp = 0;
  for (auto const& cont : *cont_frags)
  { 
    artdaq::ContainerFragment cont_frag(cont);
    numcont++;
    int numfrag = 0;
    for (size_t fragi = 0; fragi < cont_frag.block_count(); ++fragi)
    {
      numfrag++;
      artdaq::Fragment frag = *cont_frag[fragi];
      sbndaq::TDCTimestampFragment tdc_frag(frag); 
      TDCTimeStamp = tdc_frag.getTDCTimestamp()->timestamp_ns()/1e9;
    } //End of loop over the number of fragments per container
  } //End of loop over the number of containers
  return TDCTimeStamp;
}

SBNDBNBRetriever::TriggerInfo_t SBNDBNBRetriever::extractTriggerInfo(art::Event const& e) const {
  // Using TDC for current event, but PTB for previous event
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);

  art::InputTag TDC_itag("daq", "ContainerTDCTIMESTAMP");
  auto TDC_cont_frags = e.getHandle<artdaq::Fragments>(TDC_itag);

  TriggerInfo_t triggerInfo;
  double PTBTimeStamp = extractPTBTimeStamp(PTB_cont_frags);
  double TDCTimeStamp = extractTDCTimeStamp(TDC_cont_frags);

  std::cout << "PTBTimeStamp: " << PTBTimeStamp << std::endl;
  std::cout << "TDCTimeStamp: " << TDCTimeStamp << std::endl;

  triggerInfo.t_current_event = TDCTimeStamp;
  triggerInfo.t_previous_event = TDCTimeStamp-10;
  triggerInfo.number_of_gates_since_previous_event = 1;

  std::cout << std::setprecision(19) << "triggerInfo.t_current_event: " << triggerInfo.t_current_event << std::endl;
  std::cout << std::setprecision(19) << "triggerInfo.t_previous_event: " << triggerInfo.t_previous_event << std::endl;

  return triggerInfo;
}

SBNDBNBRetriever::MWRdata_t SBNDBNBRetriever::extractSpillTimes(TriggerInfo_t const& triggerInfo) const {
  
  // These lines get everything primed within the IFBeamDB
  //   They seem redundant but they are needed
  // try{auto cur_vec_temp = bfp->GetNamedVector((triggerInfo.t_previous_event)-fTimePad,"E:THCURR");} catch (WebAPIException &we) {std::cout << "caught 1" << std::endl;}      
  // try{auto cur_vec_temp = bfp->GetNamedVector((triggerInfo.t_current_event)+fTimePad,"E:THCURR");} catch (WebAPIException &we) {std::cout << "caught 1" << std::endl;}     
  try{bfp->FillCache((triggerInfo.t_current_event)+fTimePad);} catch (WebAPIException &we) {}     
  try{bfp->FillCache((triggerInfo.t_previous_event)-fTimePad);} catch (WebAPIException &we) {}      
  //std::cout << std::setprecision(19) << "bfp->GetCacheStartTime()" << bfp->GetCacheStartTime() << std::endl; 
  //std::cout << std::setprecision(19) << "bfp->GetCacheEndTime()" << bfp->GetCacheEndTime()  << std::endl; 
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
  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR Device Blocks Found : " << vars.size() << std::endl;
  // Tracking the time from the IFBeamDB
  double time_for_mwr;    
  
  // this is an iterator to track which of the 
  // three devices we will be working with
  int dev = 0;
  
  // The MWR devices are annoying and have confusing buffer
  // what we'll do is sort through all of them first and then 
  // match them to the closest spills in time
  // 

  //  int t_steps = int(((triggerInfo.t_previous_event - fTimePad) - (triggerInfo.t_current_event + fTimePad))/0.5)+25;
  int t_steps = int(((triggerInfo.t_current_event + fTimePad) - (triggerInfo.t_previous_event - fTimePad - 20.))/0.5)+25;
  mf::LogDebug("SBNDBNBRetriever") << " t_steps " << t_steps << std::endl;

  for(int t = 0; t < t_steps; t++){//Iterate through time increments
    for (std::string const& var : vars) {// Iterate through the devices
      
      //Make sure we have a device
      if(var.empty()){ 
	//mf::LogDebug("SBNDBNBRetriever") << " NO MWR DEVICES?!" << std::endl;
	continue;
      }
      /// Check the device name and interate the double-vector index
      if(var.find("M875BB") != std::string::npos ) dev = 0;
      else if(var.find("M876BB") != std::string::npos ) dev = 1;
      else if(var.find("MMBTBB") != std::string::npos ) dev = 2;
      else{
	//mf::LogDebug("SBNDBNBRetriever") << " NOT matched to a MWR DEVICES?!" << var << std::endl;
	continue;}
      
      time_for_mwr = 0;
      
      try{

	//Pull the MWR data for the device
	// these data are "packed"
	std::vector<double> packed_MWR = bfp_mwr->GetNamedVector((triggerInfo.t_previous_event)-20.-fTimePad+double(0.5*t),var,&time_for_mwr);

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

  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[0] times : " << MWR_times[0].size() << std::endl;	
  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[0]s : " << unpacked_MWR[0].size() << std::endl;	
  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[1] times : " << MWR_times[1].size() << std::endl;	
  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[1]s : " << unpacked_MWR[1].size() << std::endl;	
  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[2] times : " << MWR_times[2].size() << std::endl;	
  mf::LogDebug("SBNDBNBRetriever") << " Number of MWR[2]s : " << unpacked_MWR[2].size() << std::endl;	
  
  return { std::move(MWR_times), std::move(unpacked_MWR) };
}

int SBNDBNBRetriever::matchMultiWireData(
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

  mf::LogDebug("SBNDBNBRetriever") << "matchMultiWireData:: Number of time spills : " << times_temps.size() << std::endl;

  // We'll keep track of how many of these spills match to our 
  // DAQ trigger times
  int spill_count = 0;
  int spills_removed = 0;
  std::vector<int> matched_MWR;
  matched_MWR.resize(3);
  
  std::cout << std::setprecision(19) << "t_previous_event: " << triggerInfo.t_previous_event << std::endl;
  std::cout << std::setprecision(19) << "t_current_event: " << triggerInfo.t_current_event << std::endl;
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

  //  mf::LogDebug("SBNDBNBRetriever") << "Total number of Times we're going to test: " << times_temps.size() <<  std::endl;
  // mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19) << "Upper Limit : " << (triggerInfo.t_current_event)+fTimePad <<  std::endl;
  // mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19) << "Lower Limit : " << (triggerInfo.t_previous_event)+fTimePad <<  std::endl;
  
  // Iterating through each of the beamline times
  for (size_t i = 0; i < times_temps.size(); i++) {
    
    // Only continue if these times are matched to our DAQ time
    // mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19) << "Time # : " <<  i << std::endl;

    if(!isFirstEventInRun){//We already addressed the "first event" above
      if(times_temps[i] > (triggerInfo.t_current_event)+fTimePad){
	//mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19) << "Removed!  : " << times_temps[i] << std::endl;
	spills_removed++; 
	continue;} 
      if(times_temps[i] <= (triggerInfo.t_previous_event)+fTimePad){
	spills_removed++; 
	//mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19) << "Removed!  : " << times_temps[i] << std::endl;
	continue;}
    }

    //check if this spill is is minbias   
    /*
      40 ms was selected to be close to but outside the 66 ms 
      time of the next spill (when the beam is running at 15 Hz) 
      DocDB 33155 provides documentation of this
    */

    // mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19) << "matchMultiWireData:: trigger type : " << get_trigger_type_matching_gate(db, callback_trigger_type, run_number, times_temps[i]*1.e9-triggerInfo.WR_to_Spill_conversion+3.6e7, 40.) << " times : spill " << times_temps[i]*1.e9 << " - " << triggerInfo.WR_to_Spill_conversion << " + " << 3.6e7 <<  std::endl;
    
    // if(get_trigger_type_matching_gate(db, callback_trigger_type, run_number, times_temps[i]*1.e9-triggerInfo.WR_to_Spill_conversion+3.6e7, 40.) == 1){
    //      mf::LogDebug("SBNDBNBRetriever") << std::setprecision(19)  << "matchMultiWireData:: Skipped a MinBias gate at : " << times_temps[i]*1000. << std::endl;

    //  continue;
    //}
      
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
    
    sbn::BNBSpillInfo spillInfo = makeBNBSpillInfo(eventID, times_temps[i], MWRdata, matched_MWR);

    beamInfos.push_back(std::move(spillInfo));

    // We do not write these to the art::Events because 
    // we can filter events but want to keep all the POT 
    // information, so we'll write it to the SubRun
    
  }//end iteration over beam device times
  
  //  mf::LogDebug("SBNDBNBRetriever") << "matchMultiWireData:: Total spills counted:  " << spill_count << "   Total spills removed : " << spills_removed <<  std::endl;

  return spill_count;
}

sbn::BNBSpillInfo SBNDBNBRetriever::makeBNBSpillInfo
  (art::EventID const& eventID, double time, MWRdata_t const& MWRdata, std::vector<int> const& matched_MWR) const
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

  try{bfp->GetNamedData(time, "E:TOR860@",&TOR860,&TOR860_time);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:TOR875",&TOR875);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:LM875A",&LM875A);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:LM875B",&LM875B);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:LM875C",&LM875C);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:HP875",&HP875);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:VP875",&VP875);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:HPTG1",&HPTG1);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:VPTG1",&VPTG1);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:HPTG2",&HPTG2);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:VPTG2",&VPTG2);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:BTJT2",&BTJT2);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}
  try{bfp->GetNamedData(time, "E:THCURR",&THCURR);}catch (WebAPIException &we) {mf::LogDebug("SBNDBNBRetriever")<< "At time : " << time << " " << "got exception: " << we.what() << "\n";}

  std::cout << std::setprecision(19) << "spill_time: " << TOR860_time << std::endl;


  std::cout << std::setprecision(13) << "TOR860: " << TOR860 << std::endl;
  std::cout << std::setprecision(13) << "TOR875: " << TOR875 << std::endl;
  std::cout << std::setprecision(13) << "LM875A: " << LM875A << std::endl;
  std::cout << std::setprecision(13) << "LM875B: " << LM875B << std::endl;
  std::cout << std::setprecision(13) << "LM875C: " << LM875C << std::endl;
  std::cout << std::setprecision(13) << "HP875: " << HP875 << std::endl;
  std::cout << std::setprecision(13) << "VP875: " << VP875 << std::endl;
  std::cout << std::setprecision(13) << "HPTG1: " << HPTG1 << std::endl;
  std::cout << std::setprecision(13) << "VPTG1: " << VPTG1 << std::endl;
  std::cout << std::setprecision(13) << "HPTG2: " << HPTG2 << std::endl;
  std::cout << std::setprecision(13) << "VPTG2: " << VPTG2 << std::endl;
  std::cout << std::setprecision(13) << "BTJT2: " << BTJT2 << std::endl;
  std::cout << std::setprecision(13) << "THCURR: " << THCURR << std::endl;

  //crunch the times 
  unsigned long int time_closest_int = (int) TOR860_time;
  double time_closest_ns = (TOR860_time - time_closest_int)*1e9;
  
  //Store everything in our data-product
  sbn::BNBSpillInfo beamInfo;
  beamInfo.TOR860 = TOR860*1e12; //add in factor of 1e12 protons to get correct POT units
  beamInfo.TOR875 = TOR875*1e12; //add in factor of 1e12 protons to get correct POT units
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

  // std::cout << "TOR860: " << TOR860 << std::endl;

  for(auto const& MWRdata: unpacked_MWR){
    std::ignore = MWRdata;
    assert(!MWRdata.empty());
  }

  if(unpacked_MWR[0].empty()){
    beamInfo.M875BB.clear();
    beamInfo.M875BB_spill_time_diff = -999;//units in seconds
  }
  else{
    beamInfo.M875BB = unpacked_MWR[0][matched_MWR[0]];
    beamInfo.M875BB_spill_time_diff = (MWR_times[0][matched_MWR[0]] - time);
  }

 if(unpacked_MWR[1].empty()){
    beamInfo.M876BB.clear();
    beamInfo.M876BB_spill_time_diff = -999;//units in seconds
 }
 else{
   beamInfo.M876BB = unpacked_MWR[1][matched_MWR[1]];
   beamInfo.M876BB_spill_time_diff = (MWR_times[1][matched_MWR[1]] - time);
 }

 if(unpacked_MWR[2].empty()){
    beamInfo.MMBTBB.clear();
    beamInfo.MMBTBB_spill_time_diff = -999;//units in seconds
  }
 else{
   beamInfo.MMBTBB = unpacked_MWR[2][matched_MWR[2]];
   beamInfo.MMBTBB_spill_time_diff = (MWR_times[2][matched_MWR[2]] - time);
 }
  // We do not write these to the art::Events because 
  // we can filter events but want to keep all the POT 
  // information, so we'll write it to the SubRun
  
  beamInfo.event = eventID.event(); // the rest of ID is known by art::SubRun
  
  return beamInfo;
}

void SBNDBNBRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

void SBNDBNBRetriever::endSubRun(art::SubRun& sr)
{
  //mf::LogDebug("SBNDBNBRetriever")<< "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
  //mf::LogDebug("SBNDBNBRetriever")<< "Total number of Selected Spills : " << fOutbeamInfos.size() << std::endl;
  
  //auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  //std::swap(*p, fOutbeamInfos);
  
  //sr.put(std::move(p), art::subRunFragment());
  
  return;
}

DEFINE_ART_MODULE(SBNDBNBRetriever)
