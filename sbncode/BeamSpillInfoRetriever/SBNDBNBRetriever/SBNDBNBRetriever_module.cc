////////////////////////////////////////////////////////////////////////
// Class:       SBNDBNBRetriever
// Plugin Type: producer 
// File:        SBNDBNBRetriever_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbncode/BeamSpillInfoRetriever/POTTools.h"

namespace sbn {
  class SBNDBNBRetriever;
}

class sbn::SBNDBNBRetriever : public art::EDProducer {
public:
  explicit SBNDBNBRetriever(fhicl::ParameterSet const & params);
  // Required functions.
  void produce(art::Event & e) override;
  void beginSubRun(art::SubRun& sr) override;
  void endSubRun(art::SubRun& sr) override;

  // Plugins should not be copied or assigned.
  SBNDBNBRetriever(SBNDBNBRetriever const &) = delete;
  SBNDBNBRetriever(SBNDBNBRetriever &&) = delete;
  SBNDBNBRetriever & operator = (SBNDBNBRetriever const &) = delete;
  SBNDBNBRetriever & operator = (SBNDBNBRetriever &&) = delete;

private:
  // Declare member data here.
  std::vector< sbn::BNBSpillInfo > fOutbeamInfos;
  double fTimePad;
  double fBESOffset;
  std::string fDeviceUsedForTiming;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp;
  std::unique_ptr<ifbeam_ns::BeamFolder> bfp_mwr;
  sbn::MWRData mwrdata;
  art::ServiceHandle<ifbeam_ns::IFBeam> ifbeam_handle;
  static constexpr double MWRtoroidDelay = -0.035; ///< the same time point is measured _t_ by MWR and _t + MWRtoroidDelay`_ by the toroid [ms]

  TriggerInfo_t extractTriggerInfo(art::Event const& e) const;
  int matchMultiWireData(
    art::EventID const& eventID, 
    TriggerInfo_t const& triggerInfo,
    MWRdata_t const& MWRdata,
    std::vector< sbn::BNBSpillInfo >& beamInfos
    ) const;
  unsigned int TotalBeamSpills;
};

sbn::SBNDBNBRetriever::SBNDBNBRetriever(fhicl::ParameterSet const & params)
  : EDProducer{params} {
  produces< std::vector< sbn::BNBSpillInfo >, art::InSubRun >();
  fTimePad = params.get<double>("TimePadding");
  fBESOffset = params.get<double>("BESOffset");
  fDeviceUsedForTiming = params.get<std::string>("DeviceUsedForTiming");
  bfp = ifbeam_handle->getBeamFolder(params.get<std::string>("Bundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("TimeWindow")));
  bfp->set_epsilon(0.02);
  bfp_mwr = ifbeam_handle->getBeamFolder(params.get<std::string>("MultiWireBundle"), params.get<std::string>("URL"), std::stod(params.get<std::string>("MWR_TimeWindow")));
  bfp_mwr->set_epsilon(0.5);
  bfp_mwr->setValidWindow(3605);
  TotalBeamSpills = 0;
}

void sbn::SBNDBNBRetriever::produce(art::Event & e)
{

  // If this is the first event in the run, then ignore it
  // We do not currently have the ability to figure out the first
  // spill that the DAQ was sensitive to, so don't try to save any
  // spill information

  if (e.event() == 1) {
    return;
  }

  mf::LogDebug("SBNDBNBRetriever")<< "ptb_event: " << e.event() << std::endl;
  TriggerInfo_t const triggerInfo = extractTriggerInfo(e);

  if (triggerInfo.t_previous_event == 0) {
    return;
  }

  TotalBeamSpills += triggerInfo.number_of_gates_since_previous_event;
  MWRdata_t const MWRdata = extractSpillTimes(triggerInfo, bfp, bfp_mwr, fTimePad, MWRtoroidDelay, mwrdata);

  int const spill_count = matchMultiWireData(e.id(), triggerInfo, MWRdata, fOutbeamInfos);

  if(spill_count != int(triggerInfo.number_of_gates_since_previous_event))
    mf::LogDebug("SBNDBNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << " \t \t ::: WRONG!"<< std::endl;
  else
    mf::LogDebug("SBNDBNBRetriever")<< "Event Spills : " << spill_count << ", DAQ Spills : " << triggerInfo.number_of_gates_since_previous_event << std::endl;
}

sbn::TriggerInfo_t sbn::SBNDBNBRetriever::extractTriggerInfo(art::Event const& e) const {
  // Using TDC for current event, but PTB for previous event. Exception for case where no TDC.
  art::InputTag PTB_itag("daq", "ContainerPTB");
  auto PTB_cont_frags = e.getHandle<artdaq::Fragments>(PTB_itag);

  art::InputTag TDC_itag("daq", "ContainerTDCTIMESTAMP");
  auto TDC_cont_frags = e.getHandle<artdaq::Fragments>(TDC_itag);

  TriggerInfo_t triggerInfo;
  PTBInfo_t PTBInfo = extractPTBInfo(PTB_cont_frags, 2);

  if (TDC_cont_frags) {
    double TDCTimeStamp = extractTDCTimeStamp(TDC_cont_frags);
    triggerInfo.t_current_event = TDCTimeStamp - fBESOffset;
  }
  else{
    // If missing TDC, use PTB instead
    mf::LogDebug("SBNDBNBRetriever") << " Missing TDC Container Fragments!!! " << std::endl;
    triggerInfo.t_current_event = PTBInfo.currPTBTimeStamp - fBESOffset;
  }

  triggerInfo.t_previous_event = PTBInfo.prevPTBTimeStamp - fBESOffset;
  triggerInfo.number_of_gates_since_previous_event = PTBInfo.GateCounter;

  double PTBandCurrOffset = PTBInfo.currPTBTimeStamp - triggerInfo.t_current_event - fBESOffset;

  // Catch for an issue seen a few times where PTB off by a second.
  // Only need to correct prevTS because either currTS is from TDC
  // or there is no offset between currTS and PTB.
  if(abs(PTBandCurrOffset) >= 0.5){
    triggerInfo.t_previous_event-=PTBandCurrOffset;
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << "Offset between PTB and TDC, " << PTBandCurrOffset << std::endl;
    mf::LogDebug("SBNDBNBZEROBIASRetriever") << "Corrected previous event TS is " << std::setprecision(19) << triggerInfo.t_previous_event << std::endl;
  }

  return triggerInfo;
}

int sbn::SBNDBNBRetriever::matchMultiWireData(
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

  mf::LogDebug("SBNDBNBRetriever") << "matchMultiWireData:: Number of time spills : " << times_temps.size() << std::endl;

  // We'll keep track of how many of these spills match to our 
  // DAQ trigger times
  int spill_count = 0;
  int spills_removed = 0;
  std::vector<int> matched_MWR;
  matched_MWR.resize(3);
  
  // Iterating through each of the beamline times
  for (size_t i = 0; i < times_temps.size(); i++) {
    if(times_temps[i] > (triggerInfo.t_current_event)+fTimePad){
      spills_removed++; 
      continue;} 
    if(times_temps[i] <= (triggerInfo.t_previous_event)+fTimePad){
      spills_removed++; 
      continue;}

    if(BrokenClock(times_temps[i], bfp)){
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
  
  return spill_count;
}

void sbn::SBNDBNBRetriever::beginSubRun(art::SubRun& sr)
{
  return;
}

void sbn::SBNDBNBRetriever::endSubRun(art::SubRun& sr)
{
  mf::LogDebug("SBNDBNBRetriever")<< "Total number of DAQ Spills : " << TotalBeamSpills << std::endl;
  mf::LogDebug("SBNDBNBRetriever")<< "Total number of Selected Spills : " << fOutbeamInfos.size() << std::endl;

  auto p =  std::make_unique< std::vector< sbn::BNBSpillInfo > >();
  std::swap(*p, fOutbeamInfos);
  
  sr.put(std::move(p), art::subRunFragment());
  
  return;
}

DEFINE_ART_MODULE(sbn::SBNDBNBRetriever)
