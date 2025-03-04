////////////////////////////////////////////////////////////////////////
// Class:       SBNDBNBZEROBIASRetriever
// Plugin Type: producer 
// File:        SBNDBNBZEROBIASRetriever_module.cc
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
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"

#include "ifdh_art/IFBeamService/IFBeam_service.h"
#include "ifbeam_c.h"
#include "sbncode/BeamSpillInfoRetriever/MWRData.h"
#include "sbncode/BeamSpillInfoRetriever/SBNDPOTTools.h"

#include "larcorealg/CoreUtils/counter.h"

namespace sbn {
  class SBNDBNBZEROBIASRetriever;
}

class sbn::SBNDBNBZEROBIASRetriever : public art::EDProducer {
public:
  explicit SBNDBNBZEROBIASRetriever(fhicl::ParameterSet const & params);
  // Required functions.
  void produce(art::Event & e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

  // Plugins should not be copied or assigned.
  SBNDBNBZEROBIASRetriever(SBNDBNBZEROBIASRetriever const &) = delete;
  SBNDBNBZEROBIASRetriever(SBNDBNBZEROBIASRetriever &&) = delete;
  SBNDBNBZEROBIASRetriever & operator = (SBNDBNBZEROBIASRetriever const &) = delete;
  SBNDBNBZEROBIASRetriever & operator = (SBNDBNBZEROBIASRetriever &&) = delete;


private:
  // Declare member data here.
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  std::vector< sbn::BNBSpillInfo > fOutbeamInfosTotal;
  double fTimePad;
  double fBESOffset;
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

  static constexpr double MWRtoroidDelay = -0.035; ///< the same time point is measured _t_ by MWR and _t + MWRtoroidDelay`_ by the toroid [ms]

  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  MWRdata_t extractSpillTimes(TriggerInfo_t const& triggerInfo) const; 
  void matchMultiWireData(
    art::EventID const& eventID, 
    TriggerInfo_t const& triggerInfo,
    MWRdata_t const& MWRdata,
    std::vector< sbn::BNBSpillInfo >& beamInfos
    ) const;
  unsigned int TotalBeamSpills;

};

sbn::SBNDBNBZEROBIASRetriever::SBNDBNBZEROBIASRetriever(fhicl::ParameterSet const & params)
  : EDProducer{params} {
  raw_data_label = params.get<std::string>("raw_data_label", "daq");
  fInputLabel = params.get<std::string>("InputLabel");
  fDeviceUsedForTiming = params.get<std::string>("DeviceUsedForTiming");
  fTimePad = params.get<double>("TimePadding");
  fBESOffset = params.get<double>("BESOffset");
  fInputNonContainerInstance = params.get<std::string>("InputNonContainerInstance");
  fOutputInstance = params.get<std::string>("OutputInstance");
  fDebugLevel = params.get<int>("DebugLevel",0);
  bfp = ifbeam_handle->getBeamFolder(params.get<std::string>("Bundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("TimeWindow")));
  bfp->set_epsilon(0.02);
  bfp_mwr = ifbeam_handle->getBeamFolder(params.get<std::string>("MultiWireBundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("MWR_TimeWindow")));
  bfp_mwr->set_epsilon(0.5);
  bfp_mwr->setValidWindow(3605);
  TotalBeamSpills = 0;
  produces< std::vector< sbn::BNBSpillInfo >, art::InEvent >();
  produces< std::vector< sbn::BNBSpillInfo >, art::InSubRun >();
}

int eventNum =0;
int _run;
int _subrun;
int _event;

void sbn::SBNDBNBZEROBIASRetriever::produce(art::Event & e)
{

  // If this is the first event in the run, then ignore it
  // We do not currently have the ability to figure out the first
  // spill that the DAQ was sensitive to, so don't try to save any
  // spill information

  if (e.event() == 1) {
    auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
    std::swap(*p, fOutbeamInfos);
    e.put(std::move(p));
    return;
  }


  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);


  if (triggerInfo.t_previous_event == 0) {
    auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
    std::swap(*p, fOutbeamInfos);
    e.put(std::move(p));
    return;
  }

  TotalBeamSpills += triggerInfo.number_of_gates_since_previous_event;
  MWRdata_t const MWRdata = extractSpillTimes(triggerInfo);

  matchMultiWireData(e.id(), triggerInfo, MWRdata, fOutbeamInfos);
  fOutbeamInfosTotal.insert(fOutbeamInfosTotal.end(),fOutbeamInfos.begin(),fOutbeamInfos.end());
  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  std::swap(*p, fOutbeamInfos);
  e.put(std::move(p));
}

sbn::TriggerInfo_t sbn::SBNDBNBZEROBIASRetriever::extractTriggerInfo(art::Event const& e) const {
  // Using TDC for current event, but PTB for previous event
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);

  art::InputTag TDC_itag("daq", "ContainerTDCTIMESTAMP");
  auto TDC_cont_frags = e.getHandle<artdaq::Fragments>(TDC_itag);

  PTBInfo_t PTBInfo;
  TriggerInfo_t triggerInfo;
  PTBInfo = extractPTBInfo(PTB_cont_frags, 1);

  if (TDC_cont_frags) {
    double TDCTimeStamp = extractTDCTimeStamp(TDC_cont_frags);
    triggerInfo.t_current_event = TDCTimeStamp - fBESOffset;
  }
  else{
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Missing TDC Container Fragments!!!" << std::endl;
    triggerInfo.t_current_event = PTBInfo.currPTBTimeStamp - fBESOffset;
  }

  triggerInfo.t_previous_event = PTBInfo.prevPTBTimeStamp - fBESOffset;
  triggerInfo.number_of_gates_since_previous_event = PTBInfo.GateCounter;

  if(triggerInfo.t_current_event + fBESOffset - PTBInfo.currPTBTimeStamp >= 1){
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << " PTB and TDC Disagree!!! Correcting by adding 1 second." << std::endl;
    triggerInfo.t_previous_event+=1;
  }
  else if(triggerInfo.t_current_event + fBESOffset - PTBInfo.currPTBTimeStamp <= -1){
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << " PTB and TDC Disagree!!! Correcting by subtracting 1 second." << std::endl;
    triggerInfo.t_previous_event-=1;
  }

  return triggerInfo;
}

sbn::MWRdata_t sbn::SBNDBNBZEROBIASRetriever::extractSpillTimes(TriggerInfo_t const& triggerInfo) const {
  
  // These lines get everything primed within the IFBeamDB.
  try{bfp->FillCache((triggerInfo.t_current_event)+fTimePad);} catch (WebAPIException &we) {};     
  try{bfp->FillCache((triggerInfo.t_previous_event)-fTimePad);} catch (WebAPIException &we) {};      
  try{bfp_mwr->FillCache((triggerInfo.t_current_event)+fTimePad);} catch (WebAPIException &we) {};
  try{bfp_mwr->FillCache((triggerInfo.t_previous_event)-fTimePad);} catch (WebAPIException &we) {};

  // The multiwire chambers provide their
  // data in a vector format but we'll have 
  // to sort through it in std::string format
  // to correctly unpack it
  std::vector< std::vector< std::vector< int > > >  unpacked_MWR;
  std::vector< std::vector< double> > MWR_times;
  unpacked_MWR.resize(3);
  MWR_times.resize(3);
  std::string packed_data_str; 
  
  // Create a list of all the MWR devices with their different
  // memory buffer increments 
  // generally in the format: "E:<Device>.{Memory Block}"
  std::vector<std::string> vars = bfp_mwr->GetDeviceList();
  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR Device Blocks Found : " << vars.size() << std::endl;
  // Tracking the time from the IFBeamDB
  double time_for_mwr;    
  
  // this is an iterator to track which of the 
  // three devices we will be working with
  int dev = 0;
  
  // The MWR devices are annoying and have confusing buffer
  // what we'll do is sort through all of them first and then 
  // match them to the closest spills in time
  // 

  int t_steps = int(((triggerInfo.t_current_event + fTimePad) - (triggerInfo.t_previous_event - fTimePad - 20.))/0.5)+25;

  for(int t = 0; t < t_steps; t++){//Iterate through time increments
    for (std::string const& var : vars) {// Iterate through the devices
      
      //Make sure we have a device
      if(var.empty()){ 
        mf::LogDebug("SBNDBNBZEROBIASRetriever") << " NO MWR DEVICES?!" << std::endl;
	continue;
      }
      /// Check the device name and interate the double-vector index
      if(var.find("M875BB") != std::string::npos ) dev = 0;
      else if(var.find("M876BB") != std::string::npos ) dev = 1;
      else if(var.find("MMBTBB") != std::string::npos ) dev = 2;
      else{
	mf::LogDebug("SBNDBNBZEROBIASRetriever") << " NOT matched to a MWR DEVICES?!" << var << std::endl;
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

  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR[0] times : " << MWR_times[0].size() << std::endl;	
  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR[0]s : " << unpacked_MWR[0].size() << std::endl;	
  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR[1] times : " << MWR_times[1].size() << std::endl;	
  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR[1]s : " << unpacked_MWR[1].size() << std::endl;	
  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR[2] times : " << MWR_times[2].size() << std::endl;	
  mf::LogDebug("SBNDBNBZEROBIASRetriever") << " Number of MWR[2]s : " << unpacked_MWR[2].size() << std::endl;	
  
  return { std::move(MWR_times), std::move(unpacked_MWR) };
}

void sbn::SBNDBNBZEROBIASRetriever::matchMultiWireData(
  art::EventID const& eventID,
  TriggerInfo_t const& triggerInfo,
  MWRdata_t const& MWRdata,
  std::vector< sbn::BNBSpillInfo >& beamInfos
) const {
  
  auto const& [ MWR_times, unpacked_MWR ] = MWRdata; // alias
  
  //Here we will start collecting all the other beamline devices
  // First we get the times that the beamline device fired
  //  we have to pick a specific variable to use
  std::vector<double> times_temps = bfp->GetTimeList(fDeviceUsedForTiming);

  mf::LogDebug("SBNDBNBZEROBIASRetriever") << "matchMultiWireData:: Number of time spills : " << times_temps.size() << std::endl;

  // We'll keep track of how many of these spills match to our 
  // DAQ trigger times
  int spills_removed = 0;
  std::vector<int> matched_MWR;
  matched_MWR.resize(3);

  // Iterating through each of the beamline times
 
  double best_diff = 10000000000.0; 
  double diff; 
  size_t i = 0;

  for (size_t k = 0; k < times_temps.size(); k++){
    diff = (triggerInfo.t_current_event + fTimePad) - times_temps[k];// diff is greater than zero!
    if( diff > 0 and diff < best_diff){
      if(times_temps[k] > (triggerInfo.t_current_event)+fTimePad){
        spills_removed++; 
        continue;} 
      if(times_temps[k] <= (triggerInfo.t_previous_event)+fTimePad){
        spills_removed++; 
        continue;}

      if(BrokenClock(times_temps[i], bfp)){
        continue;
      }

      best_diff = diff; 
      i = k;
    }
  }

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
    
}

void sbn::SBNDBNBZEROBIASRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

void sbn::SBNDBNBZEROBIASRetriever::endSubRun(art::SubRun& sr)
{
  mf::LogDebug("SBNDBNBZEROBIASRetriever")<< "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
  mf::LogDebug("SBNDBNBZEROBIASRetriever")<< "Total number of Selected Spills : " << fOutbeamInfosTotal.size() << std::endl;

  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  std::swap(*p, fOutbeamInfosTotal);

  sr.put(std::move(p), art::subRunFragment());

  return;
}

DEFINE_ART_MODULE(sbn::SBNDBNBZEROBIASRetriever)
